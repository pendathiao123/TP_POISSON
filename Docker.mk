# Docker.mk - Contient les règles Docker

DOCKER_IMAGE_NAME = my_program_image
DOCKER_CONTAINER_NAME = my_container
DOCKER_FILE = ./docker/Dockerfile

# Règle pour construire l'image Docker
docker-build:
	docker build -t $(DOCKER_IMAGE_NAME) -f $(DOCKER_FILE) .

# Règle pour exécuter l'image Docker
docker-run:
	docker run --name $(DOCKER_CONTAINER_NAME) -d $(DOCKER_IMAGE_NAME)

# Règle pour arrêter et supprimer le conteneur
docker-stop:
	docker stop $(DOCKER_CONTAINER_NAME)
	docker rm $(DOCKER_CONTAINER_NAME)

# Règle pour supprimer l'image Docker
docker-rm:
	docker rmi $(DOCKER_IMAGE_NAME)
