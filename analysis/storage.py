"""Atomic keyed storage that retains completed and failed original samples."""
import sqlite3
import json
import zlib

class Store:
    """Persist one terminal outcome per cell and replicate under a fixed design hash."""

    def __init__(self, path, design_hash, keys):
        self.db = sqlite3.connect(path)
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.execute('PRAGMA synchronous=FULL')
        self.db.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)')
        old = self.db.execute("SELECT value FROM metadata WHERE key='design_hash'").fetchone()
        if old and old[0] != design_hash:
            raise ValueError('checkpoint design hash mismatch')
        self.db.execute('INSERT OR IGNORE INTO metadata VALUES (?,?)', ('design_hash', design_hash))
        self.db.execute('CREATE TABLE IF NOT EXISTS jobs (cell INTEGER,rep INTEGER,status TEXT NOT NULL,payload BLOB,PRIMARY KEY(cell,rep))')
        self.db.executemany('INSERT OR IGNORE INTO jobs VALUES (?,?,?,NULL)', [(c, r, 'NOT_ATTEMPTED') for c, r in keys])
        self.db.commit()
        assert set(((x[0], x[1]) for x in self.db.execute('SELECT cell,rep FROM jobs'))) == set(keys)

    def status(self, c, r):
        return self.db.execute('SELECT status FROM jobs WHERE cell=? AND rep=?', (c, r)).fetchone()[0]

    def start(self, c, r):
        state = self.status(c, r)
        if state not in ('NOT_ATTEMPTED', 'RUNNING'):
            raise ValueError('completed/failed key cannot be replaced')
        with self.db:
            self.db.execute('UPDATE jobs SET status=? WHERE cell=? AND rep=?', ('RUNNING', c, r))
        return state == 'RUNNING'

    def finish(self, c, r, status, payload):
        if self.status(c, r) != 'RUNNING':
            raise ValueError('finish requires a running key')
        blob = zlib.compress(json.dumps(payload, allow_nan=False, separators=(',', ':')).encode(), 6)
        with self.db:
            self.db.execute('UPDATE jobs SET status=?,payload=? WHERE cell=? AND rep=?', (status, blob, c, r))

    def get(self, c, r):
        row = self.db.execute('SELECT status,payload FROM jobs WHERE cell=? AND rep=?', (c, r)).fetchone()
        return (row[0], json.loads(zlib.decompress(row[1])) if row[1] is not None else None)

    def records(self):
        for c, r, s, b in self.db.execute('SELECT cell,rep,status,payload FROM jobs ORDER BY cell,rep'):
            yield (c, r, s, json.loads(zlib.decompress(b)) if b is not None else None)

    def close(self):
        self.db.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        self.db.close()
