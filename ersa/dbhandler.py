"""Database handler"""
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker, backref
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists
from ersa.ersa_LL import Estimate
from ersa.parser import SharedSegment

Base = declarative_base()


class Result(Base):
    """ Table that holds non-vector result values """
    __tablename__ = 'result'
    id = Column(Integer, primary_key=True)
    indv1 = Column(String(250), nullable=False)
    indv2 = Column(String(250), nullable=False)
    d_est = Column(Integer, nullable=True)
    rel_est = Column(String(250), nullable=True)
    n = Column(Integer, nullable=False)
    total_cM = Column(Float, nullable=False)


class Likelihood(Base):
    """ Table that holds log-likelihoods for a histogram """
    __tablename__ = 'likelihood'
    id = Column(Integer, primary_key=True)
    result_id = Column(Integer, ForeignKey('result.id'))
    result = relationship(Result, backref=backref('LLs', cascade='all', uselist=True))
    d = Column(Integer, nullable=False)
    LL = Column(Float, nullable=False)


class Segment(Base):
    """ Table that holds matched segment start and end locations """
    __tablename__ = 'segment'
    id = Column(Integer, primary_key=True)
    result_id = Column(Integer, ForeignKey('result.id'))
    result = relationship(Result, backref=backref('segments', cascade='all', uselist=True))
    chromosome = Column(Integer, nullable=False)
    bp_start = Column(Integer, nullable=False)
    bp_end = Column(Integer, nullable=False)


def _init_db(url):
    """ Create all tables in the engine. """
    engine = create_engine(url)
    Base.metadata.create_all(engine)


def _connect(url):
    """ Connect to a database and create tables if it doesn't exist """
    if not database_exists(url):
        _init_db(url)
    engine = create_engine(url)
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()
    return session


def insert(url, est, seg_list):
    assert isinstance(est, Estimate)
    assert isinstance(seg_list[0], SharedSegment)
    session = _connect(url)

    d_est = est.d if est.reject else None
    # TODO handle tuple for rel_est, e.g. (Child, Parent)
    res = Result(indv1=est.indv1, indv2=est.indv2, d_est=d_est, rel_est=est.rel_est,
                 n=len(est.s), total_cM=sum(est.s))
    session.add(res)

    for alt in est.alts:
        new_LL = Likelihood(d=alt[0], LL=alt[2])
        new_LL.result = res
        session.add(new_LL)

    for seg in seg_list:
        new_seg = Segment(chromosome=seg.chrom, bp_start=seg.bpStart, bp_end=seg.bpEnd)
        new_seg.result = res

    session.commit()

def clear_one(url, indv1, indv2):
    s = _connect(url)
    print("Query # pre delete:")
    print(s.query(Result).count())


    s.query(Result).filter(Result.indv1 == indv1).filter(Result.indv2 == indv2).delete()
    s.query(Result).filter(Result.indv1 == indv2).filter(Result.indv2 == indv1).delete()

    print("")
    print(s.query(Result).count())

    s.commit()


def clear_all(url):
    session = _connect(url)
    session.query(Result).delete()
    session.query(Likelihood).delete()
    session.query(Segment).delete()
    session.commit()

if __name__ == '__main__':
    url = 'sqlite:///ersa_results.db'
    clear_one(url, 'TestA', 'TestB')
    # clear_all(url)
