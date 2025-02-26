import sqlalchemy as sa
import sqlalchemy.orm as so
from app import app, db

## include tables for shell queries to db
from app.models import User

@app.shell_context_processor
def make_shell_context():
    return {
        'sa' :sa,
        'so': so,
        'db': db,
        'app': app,
        'User': User
    }

