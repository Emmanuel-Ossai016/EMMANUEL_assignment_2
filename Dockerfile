FROM ubuntu:latest

RUN apt update; apt upgrade -y

RUN apt install tree ranger -y

CMD tree