import re
import pandas as pd

def count_patterns(filename):
    patterns = ["Rpv1_", "Rpv1-12_", "Rpv1-12-3_", "Suscpt_"]
    
    # Define dataframes
    count_df = pd.DataFrame(columns=['ID'] + patterns)
    boolean_df = pd.DataFrame(columns=['ID'] + patterns)
    
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            ID = 'OG_' + str(i + 1)
            count_row = [ID]
            boolean_row = [ID]
            
            for pattern in patterns:
                count = len(re.findall(pattern, line))
                count_row.append(count)
                
                # If pattern count is greater than 0, append 1, else append 0
                boolean_row.append(int(count > 0))
            
            count_df.loc[i] = count_row
            boolean_df.loc[i] = boolean_row
    
    return count_df, boolean_df

# Test the function with a file
count_df, boolean_df = count_patterns('All_deNovo_bdBlP_mclI2.txt')

# To save the dataframe as csv files
count_df.to_csv('pattern_counts.csv', index=False)
boolean_df.to_csv('pattern_boolean.csv', index=False)

