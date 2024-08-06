#!/bin/sh

# Define the Docker image to use
DOCKER_IMAGE='openfoam/openfoam7-paraview56'

# Set the mount directory to the current working directory
MOUNT_DIR=$(pwd)

# Define the home directory inside the container
HOME_DIR='/home/openfoam'

# Get the current user ID and group ID to ensure proper permissions
USER_ID=$(id -u)
GROUP_ID=$(id -g)

# Print information about the script and user running it
echo "Launching $DOCKER_IMAGE"
echo "User: \"$(id -un)\" (ID $USER_ID, group ID $GROUP_ID)"

# Run the Docker container with the specified options
docker run -it --rm \
    -e DISPLAY=$DISPLAY \
    -u $USER_ID:$GROUP_ID \
    -v /tmp/.X11-unix:/tmp/.X11-unix \
    -v $MOUNT_DIR:$HOME_DIR \
    $DOCKER_IMAGE

# "-v /tmp/.X11-unix:/tmp/.X11-unix" - Mount X11 socket for GUI applications (e.g. Paraview)
