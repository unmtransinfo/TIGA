export IMAGE="unmtransinfo/tiga-shiny"
export FILE="Dockerfile.shiny"
DOCKER_BUILDKIT=1 docker build . -t $IMAGE -f $FILE --no-cache --progress=plain
docker run -p 3848:3838 $IMAGE
