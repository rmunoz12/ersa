""" Database interface and management classes """
#   Copyright (c) 2015 by
#   Richard Munoz <rmunoz@nygenome.org>
#   Jie Yuan <jyuan@nygenome.org>
#   Yaniv Erlich <yaniv@nygenome.org>
#
#   All rights reserved
#   GPL license

from sqlalchemy import Column, ForeignKey, Integer, \
    String, Float, Boolean, DateTime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.pool import StaticPool
from sqlalchemy.sql import select
from sqlalchemy import create_engine
from sqlalchemy_utils import database_exists
from datetime import datetime
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
    created_date = Column(DateTime, default=datetime.utcnow)
    deleted = Column(Boolean, nullable=False, default=False)


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


class Database:
    """
    Represents operations that can be done on a database
    holding results from ersa_LL.estimate_relation() calls.

    Best used indirectly through DbManager.

    Parameters
    ----------
    path : str
        Path to database

    shared_pool : bool
        Uses a SharedPool if true, otherwise a StaticPool
    """
    def __init__(self, path, shared_pool=False):
        if shared_pool:
            self.engine = create_engine(path, connect_args={'check_same_thread': False},
                                        poolclass=StaticPool)
        else:
            self.engine = create_engine(path)
        if not database_exists(path):
            Base.metadata.create_all(self.engine)
        Base.metadata.bind = self.engine
        self.conn = None
        self.trans = None

    def connect(self):
        """ Initiate a connection and begin a transaction """
        self.conn = self.engine.connect()
        self.trans = self.conn.begin()

    def soft_delete(self, pairs):
        """
        Soft deletes (marks a boolean flag) a list of pairs

        Parameters
        ----------
        pairs : list[str]
            List of pairs, with each individual's id separated
            by ":"
        """
        keys = []
        for p in pairs:
            indv1, indv2 = p.split(":")
            s = select([Result.__table__]). \
                where((~ Result.__table__.c.deleted) &
                      (((Result.__table__.c.indv1 == indv1) & (Result.__table__.c.indv2 == indv2)) |
                       ((Result.__table__.c.indv1 == indv2) & (Result.__table__.c.indv2 == indv1))))
            res = self.conn.execute(s)
            for row in res:
                keys.append(row[Result.__table__.c.id])
        n = 0
        if keys:
            remainder = len(keys)
            while remainder > 0:
                i = 900 if remainder > 900 else remainder

                u = Result.__table__.update(). \
                    where(Result.__table__.c.id.in_(keys[-i:])). \
                    values(deleted=True)
                u_result = self.conn.execute(u)
                n += u_result.rowcount

                remainder -= i
                del keys[-i:]
            print("marked {:,} results deleted".format(n))
        return n

    def insert(self, ests, seg_lists):
        """
        Bulk insert of records obtained from ersa_LL.estimate_relation()

        Parameters
        ----------
        ests : list[Estimate]

        seg_lists : list[list[SharedSegment]]
        """
        assert isinstance(ests[0], Estimate)
        assert isinstance(seg_lists[0][0], SharedSegment)

        for i in range(len(ests)):
            est, seg_list = ests[i], seg_lists[i]

            d_est = est.d if est.reject else None
            rel_est1 = est.rel_est[0] if est.rel_est else None
            rel_est2 = est.rel_est[1] if est.rel_est else None
            insert_result = Result.__table__.insert()
            inserted_result = self.conn.execute(insert_result, indv1=est.indv1, indv2=est.indv2,
                                                d_est=d_est, rel_est1=rel_est1, rel_est2=rel_est2,
                                                n=len(est.s), total_cM=sum(est.s))
            result_id = inserted_result.inserted_primary_key[0]

            insert_LL = Likelihood.__table__.insert()
            self.conn.execute(insert_LL,
                              [{'result_id': result_id,
                                'd': alt[0], 'LL': alt[2]}
                               for alt in est.alts])
            insert_seg = Segment.__table__.insert()
            self.conn.execute(insert_seg,
                              [{'result_id': result_id, 'chromosome': seg.chrom,
                                'bp_start': seg.bpStart, 'bp_end': seg.bpEnd}
                               for seg in seg_list])

    def delete(self):
        """
        Physically deletes any results that have previously
        been soft deleted. Corresponding likelihoods and segments
        are also removed.
        """
        s = select([Result.__table__]). \
            where(Result.__table__.c.deleted)
        res = self.conn.execute(s)
        keys = []
        for row in res:
            keys.append(row[Result.__table__.c.id])
        print("Found keys: {:,}".format(len(keys)))

        n_deleted = {'r': 0, 'l': 0, 's': 0}
        if keys:
            remainder = len(keys)
            while remainder > 0:
                i = 999 if remainder > 999 else remainder

                d = Result.__table__.delete(). \
                    where(Result.__table__.c.id.in_(keys[-i:]))
                res = self.conn.execute(d)
                n_deleted['r'] += res.rowcount

                d = Likelihood.__table__.delete(). \
                    where(Likelihood.__table__.c.result_id.in_(keys[-i:]))
                res = self.conn.execute(d)
                n_deleted['l'] += res.rowcount

                d = Segment.__table__.delete(). \
                    where(Segment.__table__.c.result_id.in_(keys[-i:]))
                res = self.conn.execute(d)
                n_deleted['s'] += res.rowcount

                remainder -= i
                del keys[-i:]
            print()
            print("{:10} Rows Deleted".format("Table"))
            print("Results \t{:,}".format(n_deleted['r']))
            print("Likelihood \t{:,}".format(n_deleted['l']))
            print("Segment \t{:,}".format(n_deleted['s']))
        return n_deleted

    def commit(self):
        """ push changes in the current transaction to the database """
        self.trans.commit()

    def rollback(self):
        """ discards changes in the current transaction """
        self.trans.rollback()

    def close(self):
        """
        Closes the connection, a new connection is needed for
        any further operations.
        """
        self.conn.close()


class DbManager:
    """
    Context manager for Database

    Parameters
    ----------
    path : str
        path to the database

    shared_pool : bool
        Uses a SharedPool if true, otherwise a StaticPool

    Example
    -------
    with DbManager('sqlite:///:memory:') as db:
        db.insert(ests, seg_lists)

    See Also
    --------
    Database
    """
    def __init__(self, path, shared_pool=False):
        self.path = path
        self.shared_pool = shared_pool

    def __enter__(self):
        self.db = Database(self.path, self.shared_pool)
        self.db.connect()
        return self.db

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self.db.rollback()
        else:
            self.db.commit()
        self.db.close()

if __name__ == '__main__':
    with DbManager('sqlite:///ersa_results.db') as db:
        db.delete()
