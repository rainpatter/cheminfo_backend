from flask import Flask, jsonify, request, session
from flask_sqlalchemy import SQLAlchemy
from flask_cors import CORS
from werkzeug.security import generate_password_hash, check_password_hash
from flask_migrate import Migrate


from ames.ames import AmesPredictor


app = Flask(__name__)
CORS(app)
app.secret_key = 'test_secret_key'


app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///test.db"
db = SQLAlchemy(app)
migrate = Migrate(app, db)


## changes to db schema - https://flask-migrate.readthedocs.io/en/latest/


class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(25), unique=True, nullable=False)
    password_hash = db.Column(db.String(80), nullable=False)

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)
    
    def check_password(self, password):
        return check_password_hash(self.password_hash, password)


@app.route('/')
def home():
    return 'App running'


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

    
# smiles sent as query strings 
# http://127.0.0.1:5000/array_test?smiles=CCCC&smiles=CCCCCCC
# allows multisearchable smiles

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

# run with python3 app.py (programmatic running)
if __name__ == '__main__':
    app.run(debug=True)
    
