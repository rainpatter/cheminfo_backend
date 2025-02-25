from flask import jsonify, request, session
from ames.ames import AmesPredictor
from app import app

## models not working
from models import User

@app.route('/')
def home():
    return 'App running'


## LOGIN
## working

@app.route('/login', methods=['POST'])
def login():
    try:
        data = request.get_json()
        username = data.get('username')
        password = data.get('password')
        user = User.query.filter_by(username=username).first()
        if user and user.check_password(password):
            session['username'] = username
            return jsonify({'message': 'login successful'}), 200
        else: 
            return jsonify({'message': 'invalid username or password'}), 401
    except Exception as e:
        return jsonify({'message': str(e)}), 500        
    

## REGISTRATION    
## working
  
@app.route('/register', methods=['POST'])
def register():
    try:
        data = request.get_json()
        username = data.get('username')
        password = data.get('password')
        user = User.query.filter_by(username=username).first()
        if user:
            return jsonify({'message': 'user already exists'}), 401
        else:
            new_user = User(username=username)
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
    
## AMES CALCULATOR
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