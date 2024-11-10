# BioMuta pipeline

## Overview
The BioMuta pipeline gathers mutation data from various sources and combines them into a single dataset under common field structure.

The sources included in BioMuta are:
- **[cBioPortal](https://www.cbioportal.org)**

BioMuta gathers mutation data for the following cancers:
- DOID:4045 / muscle cancer
- DOID:10283 / prostate cancer
- DOID:3565 / meningioma
- DOID:3277 / thymus cancer
- DOID:5041 / esophageal cancer
- DOID:263 / kidney cancer
- DOID:2394 / ovarian cancer
- DOID:175 / vascular cancer
- DOID:9256 / colorectal cancer
- DOID:4606 / bile duct cancer
- DOID:11934 / head and neck cancer
- DOID:2531 / hematologic cancer
- DOID:1319 / brain cancer
- DOID:1785 / pituitary cancer
- DOID:9253 / gastrointestinal stromal tumor
- DOID:5158 / pleural cancer
- DOID:184 / bone cancer
- DOID:1612 / breast cancer
- DOID:11239 / appendix cancer
- DOID:1781 / thyroid cancer
- DOID:2174 / ocular cancer
- DOID:0060073 / lymphatic system cancer
- DOID:10534 / stomach cancer
- DOID:8618 / oral cavity cancer
- DOID:3953 / adrenal gland cancer
- DOID:1793 / pancreatic cancer
- DOID:1192 / peripheral nervous system neoplasm
- DOID:2998 / testicular cancer
- DOID:1324 / lung cancer
- DOID:3121 / gallbladder cancer
- DOID:4159 / skin cancer
- DOID:3571 / liver cancer
- DOID:363 / uterine cancer
- DOID:3070 / malignant glioma
- DOID:4362 / cervical cancer
- DOID:11054 / urinary bladder cancer

## Features
BioMuta pipeline comprises two steps:
1. Download
Downloads mutation lists from each source.
TBA: cBioPortal fields, cBioPortal studies
2. Convert
Formats all resources to the BioMuta standard for both data and field structure.

## Usage
## Project Structure
## License
## Acknowledgements
The liftover from GRCh37 to GRCh38 was performed with the [LiftOver](https://genome-euro.ucsc.edu/cgi-bin/hgLiftOver?hgsid=346179925_bFhSTQmbua17iNFkjSMh5Lou4CSU) command line tool developed by UCSC (insert link).

## Setting config parameters
After cloning this repo, you will need to set the parameters given in pipeline/config.json.



# Deprecated documentation
## Requirements
The following must be available on your server:

* Node.js and npm
* docker


### Setting config parameters
After cloning this repo, you will need to set the parameters given in
cof/config.json. The "server" paramater can be "tst" or "prd" for
test or production server respectively. The "app_port" is the port
in the host that should map to docker container for the app.


### Creating and starting docker container for the APP

From the "app" subdirectory, run the python script given to build and start container:
  ```
  python3 create_app_container.py -s {DEP}
  docker ps --all
  ```
The last command should list docker all containers and you should see the container
you created "running_hivelab_app_{DEP}". To start this container, the best way is
to create a service file (/usr/lib/systemd/system/docker-hivelab-app-{DEP}.service),
and place the following content in it.

  ```
  [Unit]
  Description=Glyds APP Container
  Requires=docker.service
  After=docker.service

  [Service]
  Restart=always
  ExecStart=/usr/bin/docker start -a running_hivelab_app_{DEP}
  ExecStop=/usr/bin/docker stop -t 2 running_hivelab_app_{DEP}

  [Install]
  WantedBy=default.target
  ```
This will allow you to start/stop the container with the following commands, and ensure
that the container will start on server reboot.

  ```
  $ sudo systemctl daemon-reload 
  $ sudo systemctl enable docker-hivelab-app-{DEP}.service
  $ sudo systemctl start docker-hivelab-app-{DEP}.service
  $ sudo systemctl stop docker-hivelab-app-{DEP}.service
  ```


### Mapping APP and API containers to public domains
To map the APP and API containers to public domains (e.g. www.hivelab.org and api.hivelab.org),
add apache VirtualHost directives. This VirtualHost directive can be in a new f
ile (e.g. /etc/httpd/conf.d/hivelab.conf).

  ```
  <VirtualHost *:443>
    ServerName www.hivelab.org
    ProxyPass / http://127.0.0.1:{APP_PORT}/
    ProxyPassReverse / http://127.0.0.1:{APP_PORT}/
  </VirtualHost>

  ```

where {APP_PORT} and {API_PORT} are your port for the APP and API ports 
in conf/config.json file. You need to restart apache after this changes using 
the following command:

   ```
   $ sudo apachectl restart 
   ```





