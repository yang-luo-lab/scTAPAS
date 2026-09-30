process ChunkLimits {
    label 'sctapas_env'
    input:
        tuple val(chrom), path(chunks_csv)

    output:
        path "aug_chunks_${chrom}.csv"

    script:
    """
    python3 - <<'EOF'
import pandas as pd
chunks = pd.read_csv("${chunks_csv}", header=None, names=["chrom", "start", "end"])
chunks = chunks.sort_values("start").reset_index(drop=True)
chunks["vcf_name"] = (chunks["chrom"].astype(str) + "." + chunks["start"].astype(str)
                      + "." + chunks["end"].astype(str) + ".quilt2.diploid.vcf.gz")
mids = (chunks["end"].values[:-1] + chunks["start"].values[1:]) // 2
chunks["start_bound"] = [chunks["start"].iloc[0]] + list(mids + 1)
chunks["end_bound"]   = list(mids) + [chunks["end"].iloc[-1]]
chunks[["vcf_name", "start_bound", "end_bound", "chrom"]].to_csv(
    "aug_chunks_${chrom}.csv", index=False)
EOF
    """
}
