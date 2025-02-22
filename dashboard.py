import json
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from dash import Dash, dcc, html
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output, State
from plotly.subplots import make_subplots

TARGET_PTM = {'Citrullination': {'name': 'Citrullination',
                                 'path': 'data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib.parquet',
                                 'symbol': r'UniMod:7', 'site': 'R', 'mass_shift': 0.984016}}

# Initialize the Dash app
app = Dash(__name__)

# Define the layout
app.layout = html.Div([
    dcc.Dropdown(
        id='file-dropdown',
        options=[{'label': key, 'value': key} for key in TARGET_PTM.keys()],
        value='Citrullination',
        placeholder='Select File'
    ),
    dcc.Dropdown(
        id='cancer-type-dropdown',
        options=[],
        value=None,
        multi=True,
        placeholder='Select Cancer Types'
    ),
    dcc.Graph(id='histogram-plot')
])

@app.callback(
    [Output('cancer-type-dropdown', 'options')],
    [Input('file-dropdown', 'value')]
)
def update_dropdowns(selected_ptm):
    if not selected_ptm:
        return [[]]
    ptm_info = TARGET_PTM[selected_ptm]
    df = pd.read_parquet(ptm_info['path'])
    if df is None:
        return [[]]
    global total_intensity
    metadata = pd.read_csv("data/metadata.csv")
    dat = df.merge(metadata, on='ID', how='left')
    dat['has_target_PTM'] = dat['Modified.Sequence'].str.contains(ptm_info['symbol'])
    dat['has_target_site'] = dat['Stripped.Sequence'].str.contains(ptm_info['site'])
    dat['num_PTM'] = dat['Modified.Sequence'].str.count(ptm_info['symbol'])
    total_intensity = dat.groupby(['ID', 'assay', 'has_target_PTM']).agg(total_intensity=('Precursor.Quantity', 'sum')).reset_index()
    total_intensity['total_intensity'] = np.log10(total_intensity['total_intensity'])
    total_intensity = total_intensity.merge(metadata, on='ID')
    cancer_types = [{'label': ct, 'value': ct} for ct in total_intensity['Cancer Type'].dropna().unique()]
    return [cancer_types]

@app.callback(
    Output('histogram-plot', 'figure'),
    [Input('cancer-type-dropdown', 'value'),
     State('file-dropdown', 'value')]
)
def update_histogram(selected_cancer_types, selected_ptm):
    if not selected_cancer_types:
        return go.Figure()
    ptm_info = TARGET_PTM[selected_ptm]
    total_intensity_labeled = total_intensity[total_intensity['has_target_PTM']].dropna(subset=['Cancer Type'])
    total_intensity_labeled = total_intensity_labeled[total_intensity_labeled['Cancer Type'].isin(selected_cancer_types)]
    total_intensity_labeled = total_intensity_labeled[total_intensity_labeled['assay'].isin(['FT', 'IgB'])]
    unique_assays = total_intensity_labeled['assay'].unique()
    fig = make_subplots(rows=2, cols=1, subplot_titles=unique_assays)
    cmap = plt.get_cmap('viridis')
    unique_cancer_types = total_intensity_labeled['Cancer Type'].unique()
    colors = {ct: matplotlib.colors.rgb2hex(cmap(i / len(unique_cancer_types))) for i, ct in enumerate(unique_cancer_types)}
    for i, assay in enumerate(unique_assays):
        subset_assay = total_intensity_labeled[total_intensity_labeled['assay'] == assay]
        for ct in unique_cancer_types:
            subset = subset_assay[subset_assay['Cancer Type'] == ct]
            fig.add_trace(go.Histogram(
                x=subset['total_intensity'],
                name=ct,
                marker_color=colors[ct],
                opacity=0.5
            ), row=(i + 1), col=1)
    fig.update_layout(
        barmode='overlay',
        xaxis_title='Total intensity (log10)',
        yaxis_title='Frequency',
        title=f"Distribution of total intensity of modified peptide with {ptm_info['name']} on site {ptm_info['site']} by Cancer Type",
        legend_title_text='Cancer Type'
    )
    return fig

# Run the app
if __name__ == '__main__':
    app.run_server(debug=True)