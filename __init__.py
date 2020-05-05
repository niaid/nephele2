# -*- coding: utf-8 -*-

import traceback
from contextlib import contextmanager
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from nephele2.rds.rds_config import _RO_URI, _RW_URI


@contextmanager
def db_write():
    write_engine = create_engine(_RW_URI)
    sess_mkr = sessionmaker(bind=write_engine)
    session = sess_mkr()
    try:
        yield session
        session.commit()
    except:
        print('write except')
        print(str(traceback.format_exc()))
        session.rollback()
        raise
    finally:
        session.close()


@contextmanager
def db_read():
    """
    Creates a read context connection pool to the database.
    """
    read_engine = create_engine(_RO_URI)
    sess_mkr = sessionmaker(bind=read_engine)
    session = sess_mkr()
    try:
        yield session
    except:
        print('read except')
        print(str(traceback.format_exc()))
        session.rollback()
        raise
    finally:
        session.close()
        # If you dispose of the engine, you kill the pool.
        # _db.engine.dispose()
