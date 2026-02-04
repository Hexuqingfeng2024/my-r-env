FROM rocker/rstudio:latest

# 安装常用的 R 语言系统依赖，防止后续装包失败
RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    && rm -rf /var/lib/apt/lists/*
