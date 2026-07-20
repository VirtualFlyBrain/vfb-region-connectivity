import requests
import pandas as pd
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


kb = 'http://kb.virtualflybrain.org/db/neo4j/tx/commit'
kbw = 'http://kbw.virtualflybrain.org:7474/db/neo4j/tx/commit'
pdb = 'http://pdb.virtualflybrain.org/db/neo4j/tx/commit'
pdb_dev = 'http://pdb-dev.virtualflybrain.org/db/neo4j/tx/commit'
auth = ("neo4j", "vfb")


def _make_session(retries=5, backoff_factor=1):
    """Session that retries transient connection/server errors.

    backoff_factor=1 waits 0, 2, 4, 8, 16s between retries. allowed_methods
    includes POST (not retried by default) so dropped connections to the
    neo4j endpoint are retried rather than failing immediately."""
    retry = Retry(
        total=retries,
        connect=retries,
        read=retries,
        status=retries,
        backoff_factor=backoff_factor,
        status_forcelist=(500, 502, 503, 504),
        allowed_methods=frozenset(['POST', 'GET']),
        raise_on_status=False,
    )
    session = requests.Session()
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session


def query_neo4j(query, url=kb, auth=auth, verbose=False, timeout=120):
    """Runs cypher query and returns results as a Dataframe.
    No results gives empty dataframe, error returns False.
    Retries transient connection failures before giving up."""

    json_query = {"statements": [{"statement": query}]}

    if verbose:
        print(f'Running query: {query}')
    try:
        with _make_session() as session:
            response = session.post(
                url,
                json=json_query,
                auth=auth,
                headers={'Content-Type': 'application/json'},
                timeout=timeout,
                )
        response_json = response.json()
        if verbose:
            print("Response:", response_json)

    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return False
    
    flat_json = pd.json_normalize(response_json['results'][0]['data'])
    if flat_json.empty:
        return pd.DataFrame(columns=response_json['results'][0]['columns'])
    else:
        response_dataframe = pd.DataFrame.from_records(
            data=flat_json['row'], 
            columns=response_json['results'][0]['columns']
            )
        return response_dataframe


# FOR TESTING

if __name__ == '__main__':
    query = ("""
    MATCH (n:Class)
    WHERE n.short_form IN ['FBbt_00017015']
    RETURN n.short_form, n.label
    """)
    
    results = query_neo4j(query=query, verbose=True)
    print(results)
