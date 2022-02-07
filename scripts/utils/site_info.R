# Provide info on S-binding sites

library(tibble)

# ---- Basic info ---------------------------------------------------------------------------------
# Binding sites taken from human ACE2 genbank entry (NP_001358344.1)
S_BINDING_SITES <- tribble(
  ~start_pos, ~stop_pos, ~name,
  30,         41,        "ECO:0000269|PubMed:15791205 1",
  82,         84,        "ECO:0000269|PubMed:15791205 2",
  353,        357,       "ECO:0000269|PubMed:15791205 3"
)

ALL_S_BINDING_INDS <- mapply(seq,
                             from = S_BINDING_SITES$start_pos, 
                             to = S_BINDING_SITES$stop_pos) %>% 
  unlist()
