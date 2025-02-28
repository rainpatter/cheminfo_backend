from flask import jsonify, request, session
from flask_login import current_user, login_user
import sqlalchemy as sa
from app import app
from app import db

from ames.ames import AmesPredictor
from exposure_models.exposure_calculator_worker import calculate_all
from chembl_api.api_functions import retrieve_mol_data_from_smiles
from utils.mol_handlers import smiles_from_query_string
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
# investigate Flask-WTF for form submission - CSRF attacks
# smiles MUST BE ENTERED WITH %3 instead of = from client side

@app.route("/ames/post_smiles", methods=["POST"])
def handle_ames():
    if request.method == 'POST':
        try:
            data = request.get_json()
            smiles = data.get('smiles')
            cleaned_smiles = smiles_from_query_string(smiles) 
            ames_packet = AmesPredictor(cleaned_smiles)
            return jsonify(ames_packet), 200
        except Exception as e:
            return jsonify({'message': str(e)}), 400


# Working
# FOR POSTMAN:
# {
#             "substance_name": "ethanol",
#             "cas_number": "1111-22-3",
#             "mol_weight": 46.079,
#             "long_term_inhalation": 10,
#             "long_term_dermal": 10,
#             "short_term_inhalation": 10,
#             "local_dermal": 10,
#             "proc": "PROC7",
#             "ind_prof": "ind",
#             "phys_state": "liquid",
#             "fugacity": "very low",
#             "ventilation": "indoors - no or basic ventilation",
#             "duration": ">4hr",
#             "concentration": ">25%",
#             "lev": "no",
#             "rpe_mask": "no RPE",
#             "ppe_gloves": "no PPE",
#             "lev_dermal": "no"
#             }
# error message required for unsupported data types

@app.route("/exposure_models/worker/post_chemical", methods=["POST"])
def handle_worker_exposure_model():
    if request.method == 'POST':
        try:
            data = request.get_json()
            chemical_dict = {
                'substance_name': data.get('substance_name'),
                'cas_number': data.get('cas_number'),
                'mol_weight': data.get('mol_weight'),
                'long_term_inhalation': data.get('long_term_inhalation'),
                'long_term_dermal': data.get('long_term_dermal'),
                'short_term_inhalation': data.get('short_term_inhalation'),
                'local_dermal': data.get('local_dermal'),
                'proc': data.get('proc'),
                'ind_prof': data.get('ind_prof'),
                'phys_state': data.get('phys_state'),
                'fugacity': data.get('fugacity'),
                'ventilation': data.get('ventilation'),
                'duration': data.get('duration'),
                'concentration': data.get('concentration'),
                'lev': data.get('lev'),
                'rpe_mask': data.get('rpe_mask'),
                'ppe_gloves': data.get('ppe_gloves'),
                'lev_dermal': data.get('lev_dermal')
            }
            print(chemical_dict)
            if any(x == None for x in chemical_dict.values()):
                return jsonify({'message': 'inputs have not been completed'})
            else:
                exposure_packet = calculate_all(chemical_dict), 400
                return jsonify(exposure_packet), 200
        except Exception as e:
            return jsonify({'message': str(e)}), 500


# WORKING

@app.route('/chembl/similarity', methods=["POST"])
def handle_smiles_similarity_search():
    if request.method == 'POST':
        try:
            data = request.get_json()
            smiles = data.get('smiles')
            parsed = smiles_from_query_string(smiles)
            ames_packet = retrieve_mol_data_from_smiles(parsed)
            return jsonify(ames_packet), 200
        except Exception as e:
            return jsonify({'message': str(e)}), 400

# investigate @login_required decorator for prohibited routes
