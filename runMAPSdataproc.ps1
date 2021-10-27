hailctl dataproc start loic --region us-east1 --num-workers=2 --num-secondary-workers=10 --project loic-thibaut-dev
hailctl dataproc submit --region us-east1 loic extractcontextSFS.py
hailctl dataproc stop --region us-east1 loic
