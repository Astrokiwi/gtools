import socket


def gizmoDir():
    sname = socket.gethostname()

    if ( sname=="trillian" ):
        return "/export/1/djw/gizmos"
    elif ( sname=="srv01921" ):
        return "/srv/djw1g16/gizmos"
    
    raise Exception("Unknown server; add server and directory to gizmodatadir.py")

def movieDir():
    sname = socket.gethostname()

    if ( sname=="trillian" ):
        return "/export/1/djw/movies"
    elif ( sname=="srv01921" ):
        return "/srv/djw1g16/movies"
    
    raise Exception("Unknown server; add server and directory to gizmodatadir.py")
