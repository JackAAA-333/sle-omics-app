FROM python:3.10-slim
WORKDIR /app
COPY ["全流程非线性自动组学分析/requirements.txt", "./"]
RUN pip install --no-cache-dir -r requirements.txt
COPY . /app
EXPOSE 8501
CMD ["streamlit", "run", "全流程非线性自动组学分析/app.py", "--server.port", "8501", "--server.address", "0.0.0.0"]
