export IMAGE="unmtransinfo/tiga-shiny:shiny"
export FILE="Dockerfile.shiny"
DOCKER_BUILDKIT=1 docker build --no-cache -f $FILE -t $IMAGE .
docker run -p 3848:3838 $IMAGE
