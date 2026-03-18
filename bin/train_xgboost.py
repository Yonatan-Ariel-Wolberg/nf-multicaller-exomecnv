import pandas as pd
import xgboost as xgb
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
import numpy as np

def prepare_training_data(raw_df, num_callers=7):
    caller_flags = raw_df['supp_vec'].apply(lambda x: pd.Series(list(str(x))), convert_dtype=True)
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
