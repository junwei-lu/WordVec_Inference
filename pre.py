import numpy as np
import pandas as pd
import random
data = pd.read_parquet("rawdata.parquet")
# compute the number of patients
n = len(set(data['PatientOrder']))
# number of average time length (days) of a patient visiting the hospital
pe = random.sample(data['PatientOrder'].unique().tolist(), 1000)
result_df = data.groupby('PatientOrder')['NumDays'].agg(['max', 'min']).reset_index()
result_df['Max_Min_NumDays'] = result_df['max'] - result_df['min'] + 1
len_day = result_df['Max_Min_NumDays'].mean()
# average number of words per patient per day
result_df2 = data['PatientOrder'].value_counts().reset_index()
result_df2.columns = ['PatientOrder', 'NumOccurrences']
result_df = pd.merge(result_df[['PatientOrder','Max_Min_NumDays']], result_df2, on='PatientOrder', how='left')
result_df['num_words'] = result_df['NumOccurrences']/result_df['Max_Min_NumDays']
num_words = result_df['num_words'].mean()
