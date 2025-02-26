from flask import jsonify, request, session
from flask_login import current_user, login_user
import sqlalchemy as sa
from app import app
from app import db

from ames.ames import AmesPredictor
# models not working


@app.route('/')
def home():
    return 'App running'


# LOGIN
# working
# currently logging req body as login_user user

@app.route('/login', methods=['POST'])
def login():
    from app.models import User
    try:
        data = request.get_json()
        username = data.get('username')
        password = data.get('password')
        user = User.query.filter_by(username=username).first()
        if user and user.check_password(password):
            session['username'] = username
            login_user(user, remember=data.remember_me.data)
            return jsonify({'message': 'login successful'}), 200
        else:
            return jsonify({'message': 'invalid username or password'}), 401
    except Exception as e:
        return jsonify({'message': str(e)}), 500


# REGISTRATION
# working

@app.route('/register', methods=['POST'])
def register():
    from app.models import User
    try:
        data = request.get_json()
        username = data.get('username')
        email = data.get('email')
        password = data.get('password')
        user = User.query.filter_by(username=username).first()
        email = User.query.filter_by(email=email).first()
        if user:
            return jsonify({'message': 'user already exists'}), 401
        elif email:
            return jsonify({'message': 'email address is already in use'}), 401
        else:
            new_user = User(username=username, email=email)
            new_user.set_password(password)
            db.session.add(new_user)
            db.session.commit()
            session['username'] = username
            return jsonify({'message': 'registration successful'}), 200
    except Exception as e:
        return jsonify({'message': str(e)}), 500


@app.route('/logout')
def logout():
    session.pop("username", None)
    return jsonify({'message': 'logout successful'})

# AMES CALCULATOR
# smiles sent as query strings
# http://127.0.0.1:5000/array_test?smiles=CCCC&smiles=CCCCCCC
# allows multisearchable smiles
# investigate Flask-WTF for form submission - CSRF attacks


@app.route("/ames/post_smiles", methods=["POST"])
def handle_array():
    if request.method == 'POST':
        try:
            all_smiles = request.args.getlist('smiles')
            if len(all_smiles) == 1:
                ames_packet = AmesPredictor(all_smiles[0])
                return jsonify(ames_packet), 200
            else:
                ames_packet = [AmesPredictor(x) for x in all_smiles]
                return jsonify(ames_packet), 200
        except Exception as e:
            return jsonify({'message': 'Invalid SMILES'}), 400



## investigate @login_required decorator for prohibited routes