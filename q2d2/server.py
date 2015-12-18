import tornado
import tornado.ioloop
import tornado.httpclient
import tornado.web
import webbrowser
import os
import io
import json
import time
import subprocess
import atexit
from q2d2 import data_type_to_study_filename, create_workflow

__UPLOADS__ = "uploads/"
RELATIVE = os.path.split(__file__)[0]
GB = 1024 * 1024 * 1024

class Index(tornado.web.RequestHandler):
    def get(self):
        self.render("static/index.html")

    def delete(self):
        type = self.get_query_argument('type')
        if type in data_type_to_study_filename.values():
            os.remove(type)

class Workflows(tornado.web.RequestHandler):
    def get(self):
        self.render("static/workflows.html")

def instantiator(port):
    class Instantiator(tornado.web.RequestHandler):
        def get(self, path):
            create_workflow(path)
            self.redirect('http://localhost:%d/notebooks/%s.md' % (port, path))
    return Instantiator

@tornado.web.stream_request_body
class Upload(tornado.web.RequestHandler):

    def post(self):
        self._file.close()
        self.finish(json.dumps({}))

    def prepare(self):
        self._boundary = b'--' + ''.join(self.request.headers.get('Content-Type').split('=')[1:]).encode('ascii')
        self._end = self._boundary + b'--'

        self._file = open(self.get_query_argument('type'), mode='wb')

    def data_received(self, data):
        start = data[:len(self._boundary)]
        is_boundary = self._boundary == start
        if(is_boundary):
            index = data.find(b'\r\n\r\n')
            data = data[index+4:]
        else:
            index = data.find(b'\r\n' + self._end)
            if index > 0:
                data = data[:index]

        self._file.write(data)




def start_server(port=4444):
    jport = port + 1
    application = tornado.web.Application([
            (r"/", Index),
            (r"/workflows", Workflows),
            (r'/static/(.*)', tornado.web.StaticFileHandler,
             {'path': RELATIVE + '/static/'}),
            (r'/workflow/(.*)', instantiator(jport)),
            (r"/upload", Upload),
            ], debug=True)

    application.listen(port, max_buffer_size=10 * GB)
    process = subprocess.Popen("jupyter notebook --no-browser --port %d --port-retries 0" % jport, shell=True)
    @atexit.register
    def shutdown():
        subprocess.call("ps -ef | awk '$3 == \"%d\" {print $2}' | xargs kill -15" % process.pid, shell=True)
        process.kill()

    webbrowser.open('http://localhost:%d' % port, new=1, autoraise=True)
    tornado.ioloop.IOLoop.instance().start()
