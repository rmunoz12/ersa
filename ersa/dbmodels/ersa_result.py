""" Result database table """
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from sqlalchemy import Column, Integer, \
    String, Float, Boolean, DateTime, BigInteger
from sqlalchemy.orm import relationship
from datetime import datetime
from .base import Base

class Result(Base):
    """ Table that holds result values except for segment info """
    __tablename__ = "ersa_result"
    id = Column(Integer, primary_key=True)
    indv1 = Column(String(250), nullable=False, index=True)
    indv2 = Column(String(250), nullable=False, index=True)
    d_est = Column(Integer, nullable=True, index=True)
    rel_est1 = Column(String(250), nullable=True)
    rel_est2 = Column(String(250), nullable=True)
    n = Column(Integer, nullable=False)
    na = Column(Integer, nullable=False)
    total_cM = Column(Float, nullable=False, index=True)
    total_bp = Column(BigInteger, nullable=False, index=True)
    LLs = Column(String, nullable=False)
    segments = relationship("Segment", backref='result', cascade="all, delete, delete-orphan")
    created_date = Column(DateTime, default=datetime.utcnow, index=True)
    deleted = Column(Boolean, nullable=False, default=False, index=True)
