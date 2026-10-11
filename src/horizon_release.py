"""Fixed-horizon releases for the existing empirical variance functional.

All inputs refer to a verified common cohort with no delayed entry or weights.
The upper endpoint equals that of a truthful feasible full release under the
same censoring model. Separate-censoring lower endpoints need not be preserved.
"""
import math
from method import from_records,from_tables
from observation_model import validate

def _horizon(tau):
    if not math.isfinite(float(tau)) or tau<0:
        raise ValueError('The fixed horizon must be finite and nonnegative.')

def _size(n):
    if int(n)!=n or n<1:
        raise ValueError('The declared cohort size must be a positive integer.')

def administratively_censor(records,tau):
    """Retain every person and event precedence at the fixed horizon."""
    _horizon(tau)
    validate(records,tau)
    return [(min(t,tau),int(e and t<=tau)) for t,e in records]

def from_horizon_records(n,p,d,tau,model='separate',check_optimality=True):
    """Complete administrative records, including exits exactly at the horizon."""
    _size(n);_horizon(tau)
    if len(p)!=n or len(d)!=n:
        raise ValueError('All declared persons must remain in each margin.')
    if any(t>tau for t,e in list(p)+list(d)):
        raise ValueError('Administrative records cannot extend past the horizon.')
    return from_records(p,d,tau,model,check_optimality)

def _tables(n,tp,td,tau,terminal=False):
    _size(n);_horizon(tau)
    for t,y,d in list(tp)+list(td):
        if not math.isfinite(float(t)) or t>tau or (terminal and t==tau):
            raise ValueError('Event rows must be within the declared prefix.')

def from_horizon_tables(n,tp,td,tau,yp,yd,model='separate',check_optimality=True):
    """Exact event-risk rows through tau and exact pre-event horizon counts.

    Events exactly at tau are retained; their row risks must agree with yp/yd.
    Censor multiplicities between boundaries follow by count subtraction.
    """
    _tables(n,tp,td,tau)
    return from_tables(n,tp,td,tau,yp,yd,model,check_optimality)

def from_terminal_tables(n,tp,td,tau,yp,yd,model='separate',check_optimality=True):
    """Event-risk rows strictly before tau, with terminal counts X>=tau.

    Event status at tau is deliberately unreleased. The last censor cells
    are optimization representatives for terminal categories, not assertions
    about the original event status of any person. Their upper bound equals
    the feasible full-release upper bound. Their lower bound refers to this
    less informative release.
    """
    _tables(n,tp,td,tau,terminal=True)
    return from_tables(n,tp,td,tau,yp,yd,model,check_optimality)

def from_terminal_records(n,head_p,head_d,tau,model='separate',check_optimality=True):
    """Exact exits strictly before tau; all remaining persons form terminal groups."""
    _size(n);_horizon(tau)
    if len(head_p)>n or len(head_d)>n:
        raise ValueError('Prefix sizes exceed the declared cohort.')
    if any(t>=tau for t,e in list(head_p)+list(head_d)):
        raise ValueError('Head records must be strictly before the horizon.')
    p=list(head_p)+[(tau,0)]*(n-len(head_p));d=list(head_d)+[(tau,0)]*(n-len(head_d))
    return from_horizon_records(n,p,d,tau,model,check_optimality)
