#!/usr/bin/env bash
cmd=(datasets summary genome taxon "${snakemake_params[taxon]}" ${snakemake_params[args]} ${snakemake_params[limit]})

for filter in ${snakemake_params[filters]}; do
    if [[ "$filter" == "reference" ]]; then
        flag="--reference"
    else
        flag="--assembly-level $filter"
    fi
    echo "Trying: ${cmd[@]} $flag" >> "${snakemake_log[0]}"
    if "${cmd[@]}" $flag > "${snakemake_output[0]}" 2>> "${snakemake_log[0]}"; then
        # Check if genome summary is not empty
        if [ -s "${snakemake_output[0]}" ] && ! grep -q '"total_count": 0' "${snakemake_output[0]}"; then
            echo "Success with $flag" >> "${snakemake_log[0]}"
            exit 0
        fi
    fi
    echo "Failed with $flag, trying next..." >> "${snakemake_log[0]}"
    sleep 1
done

echo "No genome summary found for ${snakemake_wildcards[query]}" >> "${snakemake_log[0]}"
exit 1