# ============================================================
# C2S-Scale-Gemma-2-27B Cell Type Prediction
# Colab A100, 4-bit quantization, batched inference
# 07/13/2026
# Revised version:
# - Decodes only newly generated tokens
# - Prevents prompt leakage
# - Preserves highest-ranked genes during truncation
# - Removes <ctrl100> artifacts
# - Resets invalid checkpoint predictions
# - Supports safe resume
# ============================================================

!pip install -q -U \
    transformers \
    accelerate \
    bitsandbytes \
    sentencepiece \
    tqdm \
    pandas==2.2.2

# -----------------------------
# Imports
# -----------------------------
import os
import re
import gc
import torch
import pandas as pd

from tqdm.auto import tqdm
from transformers import (
    AutoTokenizer,
    AutoModelForCausalLM,
    BitsAndBytesConfig
)

# ============================================================
# Settings
# ============================================================

MODEL_ID = "vandijklab/C2S-Scale-Gemma-2-27B"

INPUT_CSV = (
    "/content/drive/MyDrive/C2S-27B/"
    "GSE283839_ATRT_RNA_V3_published_tumor_cell_sentences_top1000.csv"
)

OUTPUT_CSV = (
    "/content/drive/MyDrive/C2S-27B/"
    "GSE283839_ATRT_RNA_V3_published_tumor_27B_cell_type_predictions.csv"
)

CHECKPOINT_CSV = (
    "/content/drive/MyDrive/C2S-27B/"
    "GSE283839_ATRT_RNA_V3_published_tumor_27B_checkpoint.csv"
)

CELL_COL = "cell_name"
SENTENCE_COL = "cell_sentence"
PREDICTION_COL = "predicted_cell_type"

ORGANISM = "Homo sapiens"
NUM_GENES = 1000

# Reduce to 4 if you encounter CUDA out-of-memory errors.
BATCH_SIZE = 8

# Cell-type labels are usually short, but allow enough room for
# labels such as "CD8-positive alpha-beta T cell".
MAX_NEW_TOKENS = 24

SAVE_EVERY_BATCHES = 10

# Maximum total prompt length.
# Keep 4096 if memory usage is a concern.
MAX_INPUT_TOKENS = 4096

TEST_MODE = False
TEST_N = 20

# When True, existing invalid predictions in the checkpoint
# will be cleared and rerun.
RESET_INVALID_PREDICTIONS = True


# ============================================================
# Reproducibility and device information
# ============================================================

torch.manual_seed(42)

print("PyTorch version:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())

if torch.cuda.is_available():
    print("GPU:", torch.cuda.get_device_name(0))
    print(
        "GPU memory:",
        round(
            torch.cuda.get_device_properties(0).total_memory
            / 1024**3,
            2
        ),
        "GB"
    )


# ============================================================
# Load and validate input data
# ============================================================

df = pd.read_csv(INPUT_CSV)

required_columns = {CELL_COL, SENTENCE_COL}
missing_columns = required_columns - set(df.columns)

if missing_columns:
    raise ValueError(
        f"Input CSV is missing required columns: {missing_columns}"
    )

df[CELL_COL] = df[CELL_COL].astype(str)

# Keep missing cell sentences as empty strings for explicit checking.
df[SENTENCE_COL] = df[SENTENCE_COL].fillna("").astype(str)

if df[CELL_COL].duplicated().any():
    duplicated_cells = df.loc[
        df[CELL_COL].duplicated(keep=False),
        CELL_COL
    ].tolist()

    raise ValueError(
        "The cell_name column contains duplicated identifiers. "
        "Each cell must have a unique identifier.\n"
        f"Examples: {duplicated_cells[:10]}"
    )

if TEST_MODE:
    df = df.head(TEST_N).copy()

print("\nInput preview:")
display(df.head())

print(f"Total cells in input: {len(df):,}")

empty_sentence_count = (
    df[SENTENCE_COL].str.strip().eq("").sum()
)

print(
    f"Cells with empty cell sentences: "
    f"{empty_sentence_count:,}"
)


# ============================================================
# Load tokenizer
# ============================================================

print("\nLoading tokenizer...")

tokenizer = AutoTokenizer.from_pretrained(
    MODEL_ID,
    trust_remote_code=True,
    clean_up_tokenization_spaces=False
)

if tokenizer.pad_token_id is None:
    tokenizer.pad_token = tokenizer.eos_token

# Left padding is appropriate for decoder-only batched generation.
tokenizer.padding_side = "left"

print("Tokenizer loaded.")
print("Padding side:", tokenizer.padding_side)
print("Pad token:", tokenizer.pad_token)
print("Pad token ID:", tokenizer.pad_token_id)
print("EOS token ID:", tokenizer.eos_token_id)


# ============================================================
# Prompt template
# ============================================================

PROMPT_PREFIX = (
    f"The following is a list of {NUM_GENES} gene names ordered by "
    f"descending expression level in a {ORGANISM} cell. "
    "Your task is to give the cell type which this cell belongs to "
    "based on its gene expression.\n"
    "Cell sentence: "
)

PROMPT_SUFFIX = (
    ".\n"
    "The cell type corresponding to these genes is:"
)


# ============================================================
# Safe prompt construction
# ============================================================

def make_prompt(
    cell_sentence,
    max_input_tokens=MAX_INPUT_TOKENS
):
    """
    Construct a prompt that fits within max_input_tokens.

    The highest-ranked genes occur at the beginning of the cell
    sentence, so truncation is applied to the end of the gene list.
    The instruction prefix and answer suffix are always retained.
    """

    cell_sentence = str(cell_sentence).strip()

    prefix_ids = tokenizer(
        PROMPT_PREFIX,
        add_special_tokens=False
    )["input_ids"]

    suffix_ids = tokenizer(
        PROMPT_SUFFIX,
        add_special_tokens=False
    )["input_ids"]

    # Reserve a small margin for model-specific special tokens.
    special_token_margin = 8

    available_gene_tokens = (
        max_input_tokens
        - len(prefix_ids)
        - len(suffix_ids)
        - special_token_margin
    )

    if available_gene_tokens <= 0:
        raise ValueError(
            "MAX_INPUT_TOKENS is too small for the prompt template."
        )

    gene_ids = tokenizer(
        cell_sentence,
        add_special_tokens=False
    )["input_ids"]

    # Retain the beginning because genes are ranked by
    # descending expression.
    if len(gene_ids) > available_gene_tokens:
        gene_ids = gene_ids[:available_gene_tokens]

    shortened_sentence = tokenizer.decode(
        gene_ids,
        skip_special_tokens=True,
        clean_up_tokenization_spaces=False
    ).strip()

    return (
        PROMPT_PREFIX
        + shortened_sentence
        + PROMPT_SUFFIX
    )


# ============================================================
# Prediction cleanup and validation
# ============================================================

PROMPT_LEAKAGE_PATTERNS = [
    "The following is a list",
    "Your task is to give",
    "Cell sentence:",
    "The cell type corresponding to these genes is",
    "ordered by descending expression",
    "based on its gene expression"
]


def clean_prediction(prediction):
    """
    Clean model control tokens and return one concise cell-type label.
    """

    if prediction is None:
        return "UNRESOLVED"

    prediction = str(prediction)

    # Remove known C2S control-token artifacts.
    prediction = prediction.replace("<ctrl100>", " ")
    prediction = prediction.replace("</ctrl100>", " ")

    # Remove other visible control tokens if emitted.
    prediction = re.sub(
        r"</?ctrl\d+>",
        " ",
        prediction,
        flags=re.IGNORECASE
    )

    # Normalize whitespace.
    prediction = prediction.replace("\r", "\n")
    prediction = prediction.strip()

    # Keep the first non-empty generated line.
    lines = [
        line.strip()
        for line in prediction.splitlines()
        if line.strip()
    ]

    if not lines:
        return "UNRESOLVED"

    prediction = lines[0]

    # Remove common answer prefixes if the model repeats them.
    answer_prefix_patterns = [
        r"^the cell type is\s*[:\-]?\s*",
        r"^cell type\s*[:\-]\s*",
        r"^answer\s*[:\-]\s*",
        r"^prediction\s*[:\-]\s*"
    ]

    for pattern in answer_prefix_patterns:
        prediction = re.sub(
            pattern,
            "",
            prediction,
            flags=re.IGNORECASE
        ).strip()

    # Stop if the model starts adding an explanation.
    explanation_separators = [
        " because ",
        " based on ",
        " characterized by ",
        " as indicated by ",
        " due to ",
        " given the expression "
    ]

    lower_prediction = prediction.lower()

    for separator in explanation_separators:
        separator_index = lower_prediction.find(separator)

        if separator_index > 0:
            prediction = prediction[:separator_index].strip()
            lower_prediction = prediction.lower()

    # Remove surrounding quotation marks.
    prediction = prediction.strip(
        " \t\n\r\"'`"
    )

    # Remove trailing punctuation.
    prediction = prediction.rstrip(
        " .,:;!?-"
    ).strip()

    if not prediction:
        return "UNRESOLVED"

    return prediction


def is_invalid_prediction(prediction):
    """
    Identify empty, erroneous, prompt-leaking, or implausible outputs.
    """

    if prediction is None:
        return True

    prediction = str(prediction).strip()

    if not prediction:
        return True

    lower_prediction = prediction.lower()

    if lower_prediction in {
        "nan",
        "none",
        "null",
        "unresolved",
        "unknown",
        "n/a"
    }:
        return True

    if lower_prediction.startswith("error:"):
        return True

    if any(
        phrase.lower() in lower_prediction
        for phrase in PROMPT_LEAKAGE_PATTERNS
    ):
        return True

    if "<ctrl" in lower_prediction:
        return True

    # Cell-type labels should not contain an entire gene sentence.
    if len(prediction) > 150:
        return True

    if len(prediction.split()) > 15:
        return True

    return False


def normalize_existing_prediction(prediction):
    """
    Clean an existing checkpoint prediction without automatically
    converting all invalid values to UNRESOLVED.
    """

    if pd.isna(prediction):
        return ""

    prediction = clean_prediction(prediction)

    if prediction == "UNRESOLVED":
        return ""

    return prediction


# ============================================================
# Load model in 4-bit
# ============================================================

print("\nLoading model in 4-bit mode...")

bnb_config = BitsAndBytesConfig(
    load_in_4bit=True,
    bnb_4bit_compute_dtype=torch.bfloat16,
    bnb_4bit_use_double_quant=True,
    bnb_4bit_quant_type="nf4"
)

model = AutoModelForCausalLM.from_pretrained(
    MODEL_ID,
    quantization_config=bnb_config,
    device_map="auto",
    torch_dtype=torch.bfloat16,
    trust_remote_code=True,
    low_cpu_mem_usage=True
)

model.eval()

print("Model loaded.")
print("Model device:", model.device)


# ============================================================
# Batched prediction
# ============================================================

def batch_predict(cell_sentences):
    """
    Predict cell types for a batch of cell sentences.

    Critically, this function decodes only the newly generated
    tokens rather than decoding the prompt plus generation.
    """

    prompts = [
        make_prompt(sentence)
        for sentence in cell_sentences
    ]

    inputs = tokenizer(
        prompts,
        return_tensors="pt",
        padding=True,
        truncation=False,
        add_special_tokens=True
    )

    # Confirm that safe prompt construction succeeded.
    batch_input_length = inputs["input_ids"].shape[1]

    if batch_input_length > MAX_INPUT_TOKENS:
        raise ValueError(
            f"Tokenized batch length is {batch_input_length}, "
            f"which exceeds MAX_INPUT_TOKENS={MAX_INPUT_TOKENS}."
        )

    inputs = {
        key: value.to(model.device)
        for key, value in inputs.items()
    }

    padded_input_length = inputs["input_ids"].shape[1]

    with torch.inference_mode():
        outputs = model.generate(
            **inputs,
            max_new_tokens=MAX_NEW_TOKENS,
            do_sample=False,
            num_beams=1,
            use_cache=True,
            pad_token_id=tokenizer.pad_token_id,
            eos_token_id=tokenizer.eos_token_id
        )

    # Decoder-only models return:
    # [padded input tokens] + [newly generated tokens]
    generated_token_ids = outputs[:, padded_input_length:]

    decoded_predictions = tokenizer.batch_decode(
        generated_token_ids,
        skip_special_tokens=True,
        clean_up_tokenization_spaces=True
    )

    predictions = []

    for raw_prediction in decoded_predictions:
        cleaned_prediction = clean_prediction(raw_prediction)

        if is_invalid_prediction(cleaned_prediction):
            cleaned_prediction = "UNRESOLVED"

        predictions.append(cleaned_prediction)

    return predictions


# ============================================================
# Prompt-length diagnostics
# ============================================================

print("\nChecking prompt lengths...")

diagnostic_sample_size = min(200, len(df))

diagnostic_prompts = [
    make_prompt(sentence)
    for sentence in df[SENTENCE_COL]
    .head(diagnostic_sample_size)
    .tolist()
]

diagnostic_lengths = [
    len(
        tokenizer(
            prompt,
            add_special_tokens=True
        )["input_ids"]
    )
    for prompt in diagnostic_prompts
]

length_series = pd.Series(diagnostic_lengths)

print(
    f"Prompt-length diagnostic sample: "
    f"{diagnostic_sample_size}"
)

print(length_series.describe())

if length_series.max() > MAX_INPUT_TOKENS:
    raise RuntimeError(
        "At least one constructed prompt exceeds "
        "MAX_INPUT_TOKENS."
    )


# ============================================================
# Initialize or load checkpoint
# ============================================================

if os.path.exists(CHECKPOINT_CSV):
    print(f"\nLoading checkpoint:\n{CHECKPOINT_CSV}")

    checkpoint_df = pd.read_csv(CHECKPOINT_CSV)

    if CELL_COL not in checkpoint_df.columns:
        raise ValueError(
            f"Checkpoint does not contain '{CELL_COL}'."
        )

    checkpoint_df[CELL_COL] = (
        checkpoint_df[CELL_COL].astype(str)
    )

    if PREDICTION_COL not in checkpoint_df.columns:
        checkpoint_df[PREDICTION_COL] = ""

    # Use the current input dataframe as the authoritative source,
    # then merge predictions from the checkpoint by cell identifier.
    checkpoint_predictions = (
        checkpoint_df[
            [CELL_COL, PREDICTION_COL]
        ]
        .drop_duplicates(
            subset=[CELL_COL],
            keep="last"
        )
    )

    result_df = df.drop(
        columns=[PREDICTION_COL],
        errors="ignore"
    ).merge(
        checkpoint_predictions,
        on=CELL_COL,
        how="left"
    )

else:
    print("\nNo checkpoint found. Starting a new run.")

    result_df = df.copy()
    result_df[PREDICTION_COL] = ""


result_df[PREDICTION_COL] = (
    result_df[PREDICTION_COL]
    .apply(normalize_existing_prediction)
)


# ============================================================
# Reset invalid existing predictions
# ============================================================

invalid_existing_mask = result_df[
    PREDICTION_COL
].apply(is_invalid_prediction)

invalid_existing_count = int(
    invalid_existing_mask.sum()
)

if RESET_INVALID_PREDICTIONS:
    result_df.loc[
        invalid_existing_mask,
        PREDICTION_COL
    ] = ""

    print(
        f"Reset invalid or unresolved checkpoint predictions: "
        f"{invalid_existing_count:,}"
    )
else:
    print(
        f"Invalid existing predictions detected but not reset: "
        f"{invalid_existing_count:,}"
    )


# Empty cell sentences cannot be predicted.
empty_sentence_mask = (
    result_df[SENTENCE_COL]
    .fillna("")
    .astype(str)
    .str.strip()
    .eq("")
)

result_df.loc[
    empty_sentence_mask,
    PREDICTION_COL
] = "ERROR: Empty cell sentence"


# ============================================================
# Determine completed and remaining cells
# ============================================================

valid_prediction_mask = (
    ~result_df[PREDICTION_COL].apply(
        is_invalid_prediction
    )
)

done_cells = set(
    result_df.loc[
        valid_prediction_mask,
        CELL_COL
    ].astype(str)
)

remaining_mask = (
    ~result_df[CELL_COL].astype(str).isin(done_cells)
    & ~empty_sentence_mask
)

remaining_df = result_df.loc[
    remaining_mask,
    [CELL_COL, SENTENCE_COL]
].copy()

print(f"\nAlready completed cells: {len(done_cells):,}")
print(f"Remaining cells to predict: {len(remaining_df):,}")


# ============================================================
# Safe checkpoint saving
# ============================================================

def save_csv_safely(dataframe, output_path):
    """
    Write through a temporary file to reduce the chance of
    corrupting a checkpoint if execution is interrupted.
    """

    output_directory = os.path.dirname(output_path)

    if output_directory:
        os.makedirs(
            output_directory,
            exist_ok=True
        )

    temporary_path = output_path + ".tmp"

    dataframe.to_csv(
        temporary_path,
        index=False
    )

    os.replace(
        temporary_path,
        output_path
    )


# Save cleaned checkpoint before starting.
save_csv_safely(
    result_df,
    CHECKPOINT_CSV
)


# ============================================================
# Run batched inference
# ============================================================

failed_batches = []

batch_starts = range(
    0,
    len(remaining_df),
    BATCH_SIZE
)

progress_bar = tqdm(
    batch_starts,
    total=(
        len(remaining_df) + BATCH_SIZE - 1
    ) // BATCH_SIZE,
    desc="Predicting cell types"
)

for batch_num, start in enumerate(
    progress_bar,
    start=1
):
    batch_df = remaining_df.iloc[
        start:start + BATCH_SIZE
    ]

    cell_names = (
        batch_df[CELL_COL]
        .astype(str)
        .tolist()
    )

    cell_sentences = (
        batch_df[SENTENCE_COL]
        .astype(str)
        .tolist()
    )

    try:
        predictions = batch_predict(
            cell_sentences
        )

        if len(predictions) != len(cell_names):
            raise RuntimeError(
                "The number of predictions does not match "
                "the number of cells in the batch."
            )

        for cell_name, prediction in zip(
            cell_names,
            predictions
        ):
            result_df.loc[
                result_df[CELL_COL].astype(str)
                == cell_name,
                PREDICTION_COL
            ] = prediction

    except torch.cuda.OutOfMemoryError as error:
        print(
            "\nCUDA out-of-memory error."
        )
        print(
            "Reduce BATCH_SIZE from "
            f"{BATCH_SIZE} to 4 or 2."
        )

        save_csv_safely(
            result_df,
            CHECKPOINT_CSV
        )

        gc.collect()
        torch.cuda.empty_cache()

        raise error

    except Exception as error:
        error_message = (
            f"{type(error).__name__}: {error}"
        )

        print(
            f"\nBatch {batch_num} failed: "
            f"{error_message}"
        )

        failed_batches.append({
            "batch_number": batch_num,
            "cell_names": cell_names,
            "error": error_message
        })

        for cell_name in cell_names:
            result_df.loc[
                result_df[CELL_COL].astype(str)
                == cell_name,
                PREDICTION_COL
            ] = f"ERROR: {error_message}"

    if (
        batch_num % SAVE_EVERY_BATCHES == 0
        or batch_num
        == (
            len(remaining_df)
            + BATCH_SIZE
            - 1
        ) // BATCH_SIZE
    ):
        save_csv_safely(
            result_df,
            CHECKPOINT_CSV
        )

        completed_now = (
            ~result_df[PREDICTION_COL]
            .apply(is_invalid_prediction)
        ).sum()

        progress_bar.set_postfix(
            completed=f"{completed_now:,}",
            failed_batches=len(failed_batches)
        )

        gc.collect()

        if torch.cuda.is_available():
            torch.cuda.empty_cache()


# ============================================================
# Final cleanup
# ============================================================

result_df[PREDICTION_COL] = (
    result_df[PREDICTION_COL]
    .apply(normalize_existing_prediction)
)

# Mark any remaining blanks explicitly.
blank_final_mask = (
    result_df[PREDICTION_COL]
    .astype(str)
    .str.strip()
    .eq("")
)

result_df.loc[
    blank_final_mask,
    PREDICTION_COL
] = "UNRESOLVED"


# ============================================================
# Final validation summary
# ============================================================

prompt_leakage_mask = result_df[
    PREDICTION_COL
].astype(str).apply(
    lambda value: any(
        phrase.lower() in value.lower()
        for phrase in PROMPT_LEAKAGE_PATTERNS
    )
)

control_token_mask = result_df[
    PREDICTION_COL
].astype(str).str.contains(
    r"</?ctrl\d+>",
    case=False,
    regex=True,
    na=False
)

error_mask = result_df[
    PREDICTION_COL
].astype(str).str.startswith(
    "ERROR:",
    na=False
)

unresolved_mask = result_df[
    PREDICTION_COL
].astype(str).eq(
    "UNRESOLVED"
)

valid_final_mask = ~(
    prompt_leakage_mask
    | control_token_mask
    | error_mask
    | unresolved_mask
)

print("\nFinal validation:")
print(f"Total cells: {len(result_df):,}")
print(
    f"Valid predictions: "
    f"{valid_final_mask.sum():,}"
)
print(
    f"Unresolved predictions: "
    f"{unresolved_mask.sum():,}"
)
print(
    f"Error predictions: "
    f"{error_mask.sum():,}"
)
print(
    f"Prompt-leakage predictions: "
    f"{prompt_leakage_mask.sum():,}"
)
print(
    f"Control-token artifacts: "
    f"{control_token_mask.sum():,}"
)


# ============================================================
# Final save
# ============================================================

save_csv_safely(
    result_df,
    OUTPUT_CSV
)

save_csv_safely(
    result_df,
    CHECKPOINT_CSV
)

print("\nDone.")
print(f"Final output saved to:\n{OUTPUT_CSV}")
print(f"Checkpoint saved to:\n{CHECKPOINT_CSV}")

print("\nPrediction preview:")
display(
    result_df[
        [CELL_COL, PREDICTION_COL]
    ].head(30)
)


# ============================================================
# Optional: display the most common predictions
# ============================================================

print("\nMost common predicted cell types:")

display(
    result_df.loc[
        valid_final_mask,
        PREDICTION_COL
    ]
    .value_counts()
    .head(30)
    .rename_axis("predicted_cell_type")
    .reset_index(name="cell_count")
)
