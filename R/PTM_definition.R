# Define PTM target information
TARGET_PTM <- list(
  Citrullination_test = list(name = "Citrullination",
                             path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib_sample.parquet",
                             symbol = "Citrullination", unimod='UniMod:7', site = "R", mass_shift = 0.984016),
  Citrullination = list(name = "Citrullination",
                        path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_citrullination_report-lib.parquet",
                        symbol = "Citrullination", unimod='UniMod:7', site = "R", mass_shift = 0.984016),
  Hypusine = list(name = "Hypusine",
                  path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_hypusine_report-lib.parquet",
                  symbol = "Hypusine", unimod='UniMod:379', site = "K", mass_shift = 87.068414),
  Deoxyhypusine = list(name = "Deoxyhypusine",
                       path = "../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_deoxyhypusine_report-lib.parquet",
                       symbol = "Deoxyhypusine", unimod='UniMod:1301', site = "K", mass_shift = 71.073499),
  Trimethylation = list(name='Trimethylation',
                        path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann37_report-lib.parquet',
                        symbol = 'Trimethyl', unimod='UniMod:37', site = 'K|R', mass_shift = 42.046950),
  
  #Dimethyl: UniMod:36,28.031300,KR
  Dimethylation = list(name='Dimethylation',
                       path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann36dimethyl_report-lib.parquet',
                       symbol = 'dimethylation', unimod='UniMod:36', site = 'K|R', mass_shift = 28.031300),
  # Acetylation:   UniMod:1,42.010565,K
  Acetylation = list(name='Acetylation',
                     path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann1acetyl_report-lib.parquet',
                     symbol = 'acetylation', unimod='UniMod:1', site = 'K', mass_shift = 42.010565),
  
  # Phospho: UniMod:21,79.966331,S/T/Y
  Phospho_S = list(name='Phospho_S',
                 path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann21phosph_S_report-lib.parquet',
                 symbol = 'phospho_S', unimod='UniMod:21', site = 'S', mass_shift = 79.966331),
  Phospho_T = list(name='Phospho_T',
                   path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann21phosph_T_report-lib.parquet',
                   symbol = 'phospho_T', unimod='UniMod:21', site = 'T', mass_shift = 79.966331),
  Phospho_Y = list(name='Phospho_Y', 
                   path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann21phosph_Y_report-lib.parquet',
                   symbol = 'phospho_Y', unimod='UniMod:21', site = 'Y', mass_shift = 79.966331)
  
  
)