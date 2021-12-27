#!/bin/sh

if [ "$DATABASE" = "postgres" ]
then
    echo "Waiting for postgres..."

    sleep 15s

    echo "PostgreSQL started"
fi

#python manage.py flush --no-input
python manage.py collectstatic

exec "$@"
