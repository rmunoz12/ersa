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
from sqlalchemy.orm import relationship, sessionmaker
from sqlalchemy.pool import StaticPool
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
    rel_est1 = Column(String(250), nullable=True)
    rel_est2 = Column(String(250), nullable=True)
    n = Column(Integer, nullable=False)
    total_cM = Column(Float, nullable=False)
    LLs = relationship("Likelihood", backref='result', cascade="all, delete, delete-orphan")
    segments = relationship("Segment", backref='result', cascade="all, delete, delete-orphan")


class Likelihood(Base):
    """ Table that holds log-likelihoods for a histogram """
    __tablename__ = 'likelihood'
    id = Column(Integer, primary_key=True)
    result_id = Column(Integer, ForeignKey('result.id'))
    d = Column(Integer, nullable=False)
    LL = Column(Float, nullable=False)


class Segment(Base):
    """ Table that holds matched segment start and end locations """
    __tablename__ = 'segment'
    id = Column(Integer, primary_key=True)
    result_id = Column(Integer, ForeignKey('result.id'))
    chromosome = Column(Integer, nullable=False)
    bp_start = Column(Integer, nullable=False)
    bp_end = Column(Integer, nullable=False)


class DBhandler:
    def __init__(self, path, shared_pool=False):
        if shared_pool:
            self.engine = create_engine(path, connect_args={'check_same_thread':False}, poolclass=StaticPool)
        else:
            self.engine = create_engine(path)
        if not database_exists(path):
            # Create all tables in engine
            # if the database doesn't exist
            Base.metadata.create_all(self.engine)
        Base.metadata.bind = self.engine
        DB_Session = sessionmaker(bind=self.engine)
        self.session = DB_Session()

    def insert(self, est, seg_list):
        """
        Insert results into a database at url.

        If the pair already exists in the database, delete
        the previous results and add in the current results.
        """
        assert isinstance(est, Estimate)
        assert isinstance(seg_list[0], SharedSegment)

        q = self.session.query(Result).\
            filter(((Result.indv1 == est.indv1) & (Result.indv2 == est.indv2)) |
                   ((Result.indv1 == est.indv2) & (Result.indv2 == est.indv1)))
        if q.count():
            self._clear_one(est.indv1, est.indv2)

        d_est = est.d if est.reject else None
        rel_est1 = est.rel_est[0] if est.rel_est else None
        rel_est2 = est.rel_est[1] if est.rel_est else None
        res = Result(indv1=est.indv1, indv2=est.indv2, d_est=d_est,
                     rel_est1=rel_est1, rel_est2=rel_est2,
                     n=len(est.s), total_cM=sum(est.s))
        for alt in est.alts:
            res.LLs.append(Likelihood(d=alt[0], LL=alt[2]))
        for seg in seg_list:
            res.segments.append(Segment(chromosome=seg.chrom, bp_start=seg.bpStart, bp_end=seg.bpEnd))
        self.session.add(res)
        self.session.commit()

    def _clear_one(self, indv1, indv2):
        s = self.session
        print("Pre-delete:")
        print("R#: {} \tL#: {} \tS#: {}"
              .format(s.query(Result).count(), s.query(Likelihood).count(), s.query(Segment).count()))

        for x, y in [(indv1, indv2), (indv2, indv1)]:
            old_res = s.query(Result).filter(Result.indv1 == x).filter(Result.indv2 == y)
            for res in old_res:
                s.delete(res)

        print("Post-delete:")
        print("R#: {} \tL#: {} \tS#: {}"
              .format(s.query(Result).count(), s.query(Likelihood).count(), s.query(Segment).count()))

        s.commit()


    def _clear_all(self):
        self.session.query(Result).delete()
        self.session.query(Likelihood).delete()
        self.session.query(Segment).delete()
        self.session.commit()

if __name__ == '__main__':
    path = 'sqlite:///ersa_results.db'
    db = DBhandler(path)
    db._clear_one('TestC', 'TestB')
    db._clear_one('TestB', 'TestA')
    # clear_all()
