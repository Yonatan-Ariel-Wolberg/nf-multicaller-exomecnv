import pandas as pd
import xgboost as xgb
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
import numpy as np

# All CNV callers supported by the pipeline.  Running all seven callers
# produces fully-populated is_{caller} and qual_norm_{caller} columns in the
# feature matrix.  Callers not provided will have NaN / 0 for their columns.
SUPPORTED_CALLERS = ('canoes', 'clamms', 'xhmm', 'gatk_gcnv', 'cnvkit', 'dragen', 'indelible')

# Minimum number of distinct callers required for training the ML classifier.
# With fewer than two callers, the concordance feature is always 1 and
# provides no discriminative signal; at least two callers must contribute to
# the merged input VCF.
MIN_CALLERS_FOR_TRAINING = 2


def validate_min_callers(caller_columns):
    """Raise ValueError when fewer than MIN_CALLERS_FOR_TRAINING callers are present.

    Parameters
    ----------
    caller_columns : list[str]
        Names of per-caller flag columns found in the training data.
        These are either ``is_{caller}`` names (feature-extraction output)
        or ``caller_{i}_flag`` names (legacy supp_vec split).

    Raises
    ------
    ValueError
        If ``len(caller_columns) < MIN_CALLERS_FOR_TRAINING``.
    """
    if len(caller_columns) < MIN_CALLERS_FOR_TRAINING:
        raise ValueError(
            f"At least {MIN_CALLERS_FOR_TRAINING} CNV callers are required for "
            f"training the classifier (concordance features are not meaningful "
            f"with fewer callers). "
            f"Found {len(caller_columns)} caller column(s): {list(caller_columns)}. "
            f"Supported callers: {list(SUPPORTED_CALLERS)}."
        )


def prepare_training_data(raw_df, num_callers=7):
    # Detect named is_* caller columns produced by feature_extraction.py.
    named_caller_cols = [c for c in raw_df.columns if c.startswith('is_')]

    if named_caller_cols:
        # Feature-extraction output: is_{caller} columns already present.
        validate_min_callers(named_caller_cols)
        return raw_df

    # Legacy path: raw SURVIVOR/Truvari output with a supp_vec column.
    validate_min_callers([f'caller_{i}_flag' for i in range(num_callers)])
    caller_flags = raw_df['supp_vec'].apply(lambda x: pd.Series(list(str(x))))
    caller_flags.columns = [f'caller_{i}_flag' for i in range(num_callers)]
    caller_flags = caller_flags.astype(int)

    feature_df = pd.concat([raw_df.drop(columns=['supp_vec']), caller_flags], axis=1)
    return feature_df

def train_validation_model(X, y):
    smote = SMOTE(random_state=42)
    X_res, y_res = smote.fit_resample(X, y)
    
    model = xgb.XGBClassifier(
        n_estimators=500,
        max_depth=4,
        learning_rate=0.05,
        scale_pos_weight=len(y[y==0]) / len(y[y==1]),
        use_label_encoder=False,
        missing=np.nan
    )
    
    return model, X_res, y_res

def cross_validate_model(model, X, y, n_splits=5):
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
    scores = cross_val_score(model, X, y, cv=skf, scoring='f1')
    
    print(f"Cross-validated F1 scores: {scores}")
    print(f"Mean F1 score: {scores.mean()}")

# Load your feature data (assumed to be in a CSV format from extraction)
# df_features = pd.read_csv('features.csv')

# Prepare features and labels
# X = prepare_training_data(df_features)
# y = df_features['truth_label']  # Assuming truth_label column exists

# Train and evaluate model
# model, X_res, y_res = train_validation_model(X, y)
# cross_validate_model(model, X_res, y_res)

# Example usage:
# Ensure to uncomment the above lines and fill in actual paths and logic as per your environment.
