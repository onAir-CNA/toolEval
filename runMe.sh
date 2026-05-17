#!/bin/bash
bash reassemble.sh

docker run -id --rm \
  -p 4321:8888 \
  -v $(pwd):/sharedFolder \
  -e JUPYTER_PASSWORD_HASH='<insert_your_password_hash_here>' \
  jupyter-r-bioinfo:final

echo "JupyterLab avviato su http://localhost:4321"
