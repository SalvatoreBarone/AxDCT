FROM ubuntu:18.04

# Update Software
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && apt-get -qq install -y apt-utils > /dev/null
RUN apt-get -qq dist-upgrade -y > /dev/null && apt-get -qq autoremove -y > /dev/null

# Copy install script
ADD ./configure_environment /opt
ADD ./scripts /opt/scripts
ADD . /root

SHELL ["/bin/bash", "-c"]

# Run iideaa installation
RUN chmod +x /opt/configure_environment
RUN /opt/configure_environment

# Cleanup
RUN echo "Cleaning up..."
RUN rm -Rf ~/llvm
RUN rm -f /opt/configure_environment
RUN rm -Rf /opt/scripts
RUN rm -Rf /opt/ParadisEO-2.0/build
RUN rm -Rf /opt/opencv
RUN rm -Rf /opt/cnl
RUN rm -Rf /opt/cmake-3.14.0-rc4

# Install zsh and oh-my-zsh
RUN echo "Installing a fancy shell..."
RUN ["apt-get", "-qq", "install", "-y", "zsh", "nano"]
RUN wget --quiet https://github.com/robbyrussell/oh-my-zsh/raw/master/tools/install.sh -O - | zsh || true
RUN sed -i "s/git/git sudo docker /g" ~/.zshrc
RUN sed -i "s/robbyrussell/af-magic/g" ~/.zshrc

# Clone AxDCT Repo (impossible while repo is private)
#RUN cd ~ && git clone https://github.com/andreaaletto/AxDCT.git

RUN echo "Configuration completed. Packaging Docker Image..."