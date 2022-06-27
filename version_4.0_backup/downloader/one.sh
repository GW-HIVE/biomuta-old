token_file="/data/projects/biomuta/downloads/v-5.0/tcga/gdc-user-token.2019-07-23T17_13_24.556Z.txt"
manifest_file="/data/projects/biomuta/downloads/v-5.0/tcga/manifest-files-grch38/ACC.txt"
out_dir="/data/projects/biomuta/downloads/v-5.0/tcga/ACC/"

/usr/bin/gdc-client download -t $token_file  -m $manifest_file -d $out_dir

