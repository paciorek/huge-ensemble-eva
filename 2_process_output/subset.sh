mth_names=("Jan" "Feb" "Mar" "Apr" "May" "Jun" "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")

for((mth=0; mth<12; mth++)); do
    Rscript -e "month <- $((mth+1))$; setwd('output_22y_40m_${mth_names[$mth]}1'); source('subset.R')"
done
