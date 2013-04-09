var connect = require('connect');
connect.createServer(
    connect.static('C:/jQuest/')
).listen(1337);