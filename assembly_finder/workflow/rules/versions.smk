rule tar_find_versions:
    output:
        tar=os.path.join(dir.versions, "tar.version"),
        find=os.path.join(dir.versions, "find.version"),
    shell:
        """
        tar --version > {output.tar}
        find --version > {output.find}
        """


rule curl_version:
    output:
        os.path.join(dir.versions, "curl.version"),
    conda:
        os.path.join(dir.env, "curl.yml")
    shell:
        """
        curl --version > {output}
        """


rule datasets_version:
    output:
        os.path.join(dir.versions, "datasets.version"),
    conda:
        os.path.join(dir.env, "datasets.yml")
    shell:
        """
        datasets version > {output}
        """


rule rsync_version:
    output:
        os.path.join(dir.versions, "rsync.version"),
    conda:
        os.path.join(dir.env, "rsync.yml")
    shell:
        """
        rsync --version > {output}
        """


rule taxonkit_csvtk_versions:
    output:
        taxonkit=os.path.join(dir.versions, "taxonkit.version"),
        csvtk=os.path.join(dir.versions, "csvtk.version"),
    conda:
        os.path.join(dir.env, "taxonkit.yml")
    shell:
        """
        taxonkit version > {output.taxonkit}
        csvtk version > {output.csvtk}
        """


rule unzip_version:
    output:
        os.path.join(dir.versions, "unzip.version"),
    conda:
        os.path.join(dir.env, "unzip.yml")
    shell:
        """
        unzip -v > {output}
        """


rule all_versions:
    input:
        rsync=os.path.join(dir.versions, "rsync.version"),
        curl=os.path.join(dir.versions, "curl.version"),
        unzip=os.path.join(dir.versions, "unzip.version"),
        tar=os.path.join(dir.versions, "tar.version"),
        find=os.path.join(dir.versions, "find.version"),
        datasets=os.path.join(dir.versions, "datasets.version"),
        taxonkit=os.path.join(dir.versions, "taxonkit.version"),
        csvtk=os.path.join(dir.versions, "csvtk.version"),
    output:
        temp(os.path.join(dir.versions, "versions.flag")),
    shell:
        """
        touch {output}
        """
