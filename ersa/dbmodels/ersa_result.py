""" Result database table """
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from sqlalchemy import Column, Integer, \
    String, Float, Boolean, DateTime
from sqlalchemy.orm import relationship
from datetime import datetime
from ersa.dbmodels.base import Base

class Result(Base):
    """ Table that holds result values except for segment info """
    __tablename__ = "ersa_result"
    id = Column(Integer, primary_key=True)
    indv1 = Column(String(250), nullable=False)
    indv2 = Column(String(250), nullable=False)
    d_est = Column(Integer, nullable=True)
    rel_est1 = Column(String(250), nullable=True)
    rel_est2 = Column(String(250), nullable=True)
    n = Column(Integer, nullable=False)
    total_cM = Column(Float, nullable=False)
    LLs = Column(String, nullable=False)
    segments = relationship("Segment", backref='result', cascade="all, delete, delete-orphan")
    created_date = Column(DateTime, default=datetime.utcnow)
    deleted = Column(Boolean, nullable=False, default=False)
