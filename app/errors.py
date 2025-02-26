from app import app, db

@app.errorhandler(401)
def unauthorized_error(error):
    return 401

@app.errorhandler(500)
def internal_error(error):
    db.session_rollback()
    return 500