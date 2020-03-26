FROM python:3

COPY requirements.txt ./

RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir smaca

ENTRYPOINT ["smaca"]

CMD ["--help"] 
