FROM rocker/rstudio:latest

# 1. 深度配置 RStudio 核心参数
# - auth-none=1: 彻底关掉登录界面
# - www-verify-user-agent=0: 关键！防止 GitHub 代理导致的浏览器重定向
# - www-frame-origin=any: 允许嵌入
# - www-root-path: 强制根目录
RUN echo "auth-none=1" >> /etc/rstudio/rserver.conf && \
    echo "www-frame-origin=any" >> /etc/rstudio/rserver.conf && \
    echo "www-verify-user-agent=0" >> /etc/rstudio/rserver.conf && \
    echo "www-allow-origin=all" >> /etc/rstudio/rserver.conf

# 2. 预设一个无需密码的 rstudio 用户环境（防止内部跳转到 auth-sign-in）
RUN echo "user-ask-password=0" >> /etc/rstudio/rsession.conf

# 3. 赋予权限
RUN chmod -R 777 /home/rstudio
