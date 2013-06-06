(function() {
  var connect;

  connect = require('connect');

  connect.createServer(connect["static"]('C:/jQuest/')).listen(1337);

}).call(this);
