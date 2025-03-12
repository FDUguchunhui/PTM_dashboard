# Define PTM target information
TARGET_PTM <- list(
  Citrullination_test = list(name = "Citrullination",
                             symbol = "citrullination_sample", unimod='UniMod:7', site = "R", mass_shift = 0.984016),
  Citrullination = list(name = "Citrullination_R",
                        symbol = "citrullination", unimod='UniMod:7', site = "R", mass_shift = 0.984016),
  Hypusine = list(name = "Hypusine_K",
                  symbol = "hypusine", unimod='UniMod:379', site = "K", mass_shift = 87.068414),
  Deoxyhypusine = list(name = "Deoxyhypusine_K",
                       symbol = "Deoxyhypusine", unimod='UniMod:1301', site = "K", mass_shift = 71.073499),
  Trimethylation = list(name='Trimethylation_KR',
                        path = '../data/RAINBOW_UNIPROT_human_revi_2024_12_19_ProteinAG_diann37_report-lib.parquet',
                        symbol = 'diann37', unimod='UniMod:37', site = 'K|R', mass_shift = 42.046950),

  #Dimethyl: UniMod:36,28.031300,KR
  Dimethylation = list(name='Dimethylation_KR',
                       symbol = 'diann36dimethyl', unimod='UniMod:36', site = 'K|R', mass_shift = 28.031300),
  # Acetylation:   UniMod:1,42.010565,K
  Acetylation = list(name='Acetylation_K',
                     symbol = 'diann1acetyl', unimod='UniMod:1', site = 'K', mass_shift = 42.010565),

  # Phospho: UniMod:21,79.966331,S/T/Y
  Phospho_S = list(name='Phospho_S',
                 symbol = 'diann21phosph_S', unimod='UniMod:21', site = 'S', mass_shift = 79.966331),
  Phospho_T = list(name='Phospho_T',
                   symbol = 'diann21phosph_T', unimod='UniMod:21', site = 'T', mass_shift = 79.966331),
  Phospho_Y = list(name='Phospho_Y',
                   symbol = 'diann21phosph_Y', unimod='UniMod:21', site = 'Y', mass_shift = 79.966331),

  glycineglycine = list(name='glycineglycine_K',
             symbol = 'diann121KGG', unimod='UniMod:121', site = 'K', mass_shift = 114.042927),

  'ADP Ribose addition' = list(name='ADP-Ribosyl_R',
             symbol = 'diann213PARyl', unimod='UniMod:213', site = 'R', mass_shift = 541.061110),

  Palmitoylation = list(name='Palmitoylation_C',
             symbol = 'diann47palmityl', unimod='UniMod:47', site = 'C', mass_shift = 238.229673),

  'Oxidation to nitro' = list(name='Nitro_T',
             symbol = 'diann354nitro', unimod='UniMod:354', site = 'T', mass_shift = 44.985078),

  nitrosylation = list(name='nitrosylation_C',
             symbol = 'diann275Snitosyl', unimod='UniMod:275', site = 'C', mass_shift = 226.029390)
)
