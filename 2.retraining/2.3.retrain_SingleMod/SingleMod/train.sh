SAMPLE="HEK293T-S1"
BASE="/data/pingyc/projects/20250601_RNA/data/HEK293T-S1/SingleMod/features"
OUT="/data/pingyc/projects/20250601_RNA/data/HEK293T-S1/SingleMod/training"
KIT="002"
GPU="0"
REP="0"
COV="20"
EPOCHS="200"
BATCH="10"
DB="/path/to/HEK293T-S1_label.tsv"   # 换成你的 m6A_database（tab分隔）

mkdir -p "$OUT"

for D in A G T; do
  for R in A G; do
    for H in A C T; do
      MOTIF="${D}${R}AC${H}"
      SEQ="${BASE}/${MOTIF}_sequence.npy"
      SIG="${BASE}/${MOTIF}_signal.npy"
      EXT="${BASE}/${MOTIF}_extra.npy"
      ODIR="${OUT}/${MOTIF}/rep${REP}"

      # 若文件不全则跳过
      if [[ ! -f "$SEQ" || ! -f "$SIG" || ! -f "$EXT" ]]; then
        echo "[SKIP] ${MOTIF} missing file(s)"
        continue
      fi

      mkdir -p "$ODIR"

      python -u /data/pingyc/projects/20250601_RNA/softwares/NaRMBench/2.retraining/2.3.retrain_SingleMod/SingleMod/SingleMod_train.py \
        -v "$KIT" \
        -s "$SAMPLE" \
        -seq "$SEQ" \
        -sig "$SIG" \
        -ext "$EXT" \
        -d "$DB" \
        -m "$MOTIF" \
        -r "$REP" \
        -c "$COV" \
        -g "$GPU" \
        -o "$ODIR" \
        > "$ODIR/training.log" 2>&1
    done
  done
done
