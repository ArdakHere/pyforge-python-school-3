from celery import Celery

celery = Celery(
    'src',
    broker='redis://redis:6379/0',
    backend='redis://redis:6379/0'
)

celery.autodiscover_tasks(['src.tasks'])

celery.conf.update(task_track_started=True)
