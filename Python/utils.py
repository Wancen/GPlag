import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler

def preprocess(houseprice, start_date, end_date, region):
    mask = (houseprice['period_end'] >= start_date) & (houseprice['period_end'] <= end_date)
    mask2 = houseprice['region_name'] == region
    filtered_df = houseprice.loc[mask&mask2]
    if filtered_df.shape[0] == 0:
        return(0)
    filtered_df['period_end'] = pd.to_datetime(filtered_df['period_end'])  # ensure 'period_end' is in datetime format

    min_date = filtered_df['period_end'].min()  # get the minimum date

    # subtract the minimum date from all dates, convert to days and add as a new column
    filtered_df['days'] = (filtered_df['period_end'] - min_date).dt.days
    filtered_df['inventory_moving_avg'] = filtered_df['inventory'].rolling(window=4).mean()
    filtered_df['price_moving_avg'] = filtered_df['median_sale_price'].rolling(window=4).mean()


    # Create a scaler object
    # scaler = MinMaxScaler()
    scaler = StandardScaler()
    # Fit and transform the 'inventory' column
    # Note: the fit_transform function expects a 2D array, so we use df[['inventory']]
    filtered_df['inventory_inverse_scaled'] = scaler.fit_transform(filtered_df[['inventory_moving_avg']])
    filtered_df['price_scaled'] = scaler.fit_transform(-filtered_df[['price_moving_avg']])
    return(filtered_df.iloc[4:])

