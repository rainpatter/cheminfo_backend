from typing import Optional
from werkzeug.security import generate_password_hash, check_password_hash
import sqlalchemy as sa
import sqlalchemy.orm as so
from app import db, login
from flask_login import UserMixin


class User(UserMixin, db.Model):
    id: so.Mapped[int] = so.mapped_column(primary_key=True)
    username: so.Mapped[str] = so.mapped_column(
        sa.String(64), index=True, unique=True)
    email: so.Mapped[str] = so.mapped_column(
        sa.String(120), index=True, unique=True)
    # can be empty or nullable
    password_hash: so.Mapped[Optional[str]] = so.mapped_column(sa.String(120))

    def set_password(self, password):
        self.password_hash = generate_password_hash(password)

    def check_password(self, password):
        return check_password_hash(self.password_hash, password)

    # printable function of class
    def __repr__(self):
        return '<User {}>'.format(self.username)

# loads session user from database


@login.user_loader
def load_user(id):
    return db.session.get(User, int(id))
