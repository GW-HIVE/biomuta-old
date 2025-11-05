# Builds the docker image 
docker build -t biomuta-file-server .

# Builds the docker container 
docker run -d -p 127.0.0.1:8087:80 -v /data/shared/repos/biomuta-old/nginx-file-server/data:/usr/share/nginx/html --name biomuta-file-server-container biomuta-file-server
