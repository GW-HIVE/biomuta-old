#delete molecular profile and sample list helper files
find /data/shared/pipelines/cbioportal/mutations -type f -name "*molecular_profiles.json" -delete
find /data/shared/pipelines/cbioportal/mutations -type f -name "*sample_lists.json" -delete
