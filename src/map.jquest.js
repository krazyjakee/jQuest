(function() {
  window.Map = (function() {
    function Map() {}

    Map.settings = {
      fileExtension: 'json',
      mapDirectory: ''
    };

    Map.mapElement = false;

    Map.mapData = false;

    Map.renderedTiles = [];

    Map.tileProperties = [];

    Map.playerTile = false;

    Map.showBlocked = false;

    Map.showPaths = false;

    Map.playerLayer = 0;

    Map.isMoving = false;

    Map.loadMap = function(targetElem, mapSource, callback) {
      this.mapElement = targetElem;
      if (typeof mapSource === 'object') {
        this.mapData = mapSource;
        this.drawMap(targetElem);
        if (callback) {
          return callback();
        }
      } else if (typeof mapSource === 'string') {
        return $.getJSON(this.settings.mapDirectory + mapSource + '.' + this.settings.fileExtension, function(json) {
          window.Map.mapData = json;
          window.Map.drawMap(targetElem);
          if (callback) {
            return callback();
          }
        });
      }
    };

    Map.unloadMap = function(callback) {
      $(this.mapElement).fadeOut('fast', function() {
        $(this.mapElement).empty();
        if (callback) {
          return callback();
        }
      });
      this.mapData = false;
      this.renderedTiles = [];
      return this.tileProperties = [];
    };

    Map.drawMap = function(targetElem) {
      var currentLayer, i, index, layer, tile, _i, _j, _len, _len1, _ref, _ref1, _results;

      _ref = this.mapData.layers;
      _results = [];
      for (index = _i = 0, _len = _ref.length; _i < _len; index = ++_i) {
        layer = _ref[index];
        if (layer.type === "tilelayer") {
          $(targetElem).append("<div class=\"layer\" id=\"layer" + index + "\" />");
          currentLayer = $("#layer" + index);
          _ref1 = layer.data;
          for (i = _j = 0, _len1 = _ref1.length; _j < _len1; i = ++_j) {
            tile = _ref1[i];
            $(currentLayer).append("<div id=\"tile" + index + "-" + i + "\" class=\"tile\" />");
            this.drawTile(index, tile, i);
          }
          $(currentLayer).css('width', this.mapData.tilewidth * this.mapData.width + 'px').css('height', this.mapData.tileheight * this.mapData.height + 'px');
          $('.tile').width(this.mapData.tilewidth).height(this.mapData.tileheight);
          $('.layer:last').find('.tile').click(this.tileClick);
          _results.push($('#maploading').remove());
        } else if (layer.type === "objectgroup") {
          if (layer.name === "Player") {
            this.playerLayer = index;
          }
          _results.push($(targetElem).append("<div class=\"layer\" id=\"layer" + index + "\" />"));
        } else {
          _results.push(void 0);
        }
      }
      return _results;
    };

    Map.drawTile = function(layer, srcTile, targetTile) {
      var heightCount, i, index, offset, properties, setData, setHeight, setWidth, target, tileset, widthCount, _i, _j, _len, _ref;

      target = $("#tile" + layer + "-" + targetTile);
      setData = false;
      _ref = this.mapData.tilesets;
      for (index = _i = 0, _len = _ref.length; _i < _len; index = ++_i) {
        tileset = _ref[index];
        if (srcTile >= tileset.firstgid) {
          setData = tileset;
        } else {
          return false;
        }
      }
      setWidth = Math.floor(setData.imagewidth / setData.tilewidth) - 1;
      setHeight = Math.floor(setData.imageheight / setData.tileheight) - 1;
      heightCount = 0;
      widthCount = 0;
      for (i = _j = 1; _j < srcTile; i = _j += 1) {
        if (widthCount < setWidth) {
          widthCount++;
        } else {
          widthCount = 0;
          heightCount++;
        }
      }
      offset = {
        x: (widthCount * setData.tilewidth) + (setData.spacing * widthCount) + setData.margin,
        y: (heightCount * setData.tileheight) + (setData.spacing * heightCount) + setData.margin
      };
      if (setData.image) {
        this.renderedTiles[srcTile] = "#tile" + layer + "-" + targetTile;
        properties = this.tileProperty(srcTile);
        if (properties) {
          this.tileProperties[targetTile] = properties;
        }
        $(target).css('background-image', "url(" + setData.image + ")");
        return $(target).css('background-position', "-" + offset.x + "px -" + offset.y + "px");
      }
    };

    Map.setFocus = function(tileId, duration) {
      var curloc, loc, maxheight, maxwidth;

      loc = this.tileIdConvert(tileId);
      curloc = $('.layer:first').position();
      maxwidth = (Map.mapData.width * 32) - $('.viewport').width();
      maxheight = (Map.mapData.height * 32) - $('.viewport').height();
      loc[0] = (loc[0] * 32) - ($('.viewport').width() / 2);
      loc[1] = (loc[1] * 32) - ($('.viewport').height() / 2);
      if (loc[1] < 0) {
        loc[1] = 0;
      }
      if (loc[0] < 0) {
        loc[0] = 0;
      }
      if (loc[1] > maxheight) {
        loc[1] = maxheight;
      }
      if (loc[0] > maxwidth) {
        loc[0] = maxwidth;
      }
      this.isMoving = true;
      $('.layer-container').stop();
      return $('.layer-container').animate({
        top: "-" + loc[1] + "px",
        left: "-" + loc[0] + "px"
      }, duration, function() {
        return this.isMoving = false;
      });
    };

    Map.tileIdConvert = function(tileInput) {
      var i, tileId, width, x, y, _i;

      if (typeof tileInput === 'object') {
        tileId = tileInput[1] * this.mapData.width - 1;
        tileId += tileInput[0];
        return tileId;
      } else {
        y = 0;
        x = 0;
        width = this.mapData.width - 1;
        for (i = _i = 0; _i < tileInput; i = _i += 1) {
          if (x < width) {
            x++;
          } else {
            x = 0;
            y++;
          }
        }
        return [x, y];
      }
    };

    Map.tileIdPosition = function(tileInput) {
      if (typeof tileInput !== 'object') {
        tileInput = this.tileIdConvert(tileInput);
      }
      return [(tileInput[0] - 1) * this.mapData.tilewidth, tileInput[1] * this.mapData.tileheight];
    };

    Map.tileProperty = function(tileId) {
      var data, index, properties, _i, _len, _ref;

      tileId = tileId - 1;
      _ref = this.mapData.tilesets;
      for (index = _i = 0, _len = _ref.length; _i < _len; index = ++_i) {
        data = _ref[index];
        if (data.tileproperties && data.tileproperties[tileId]) {
          properties = data.tileproperties[tileId];
        }
      }
      return properties;
    };

    Map.tilePropertyLogic = function(prop) {
      switch (prop.property) {
        case "door":
          this.playerTile = prop.loc.split(',');
          this.unloadMap(function() {
            window.Map.loadMap(window.Map.mapElement, prop.map);
          });
      }
    };

    Map.tileClick = function(e) {
      var index, path, paths, tileId, _i, _len;

      e.preventDefault();
      tileId = $(this).attr('id');
      tileId = tileId.substr(tileId.lastIndexOf('-') + 1);
      tileId++;
      if (e.button === 2) {
        return alert();
      } else {
        paths = window.Map.makePath(tileId);
        if (paths.length) {
          window.Character.playerMove(paths);
          if (window.Map.showPaths) {
            $("#layer" + Map.playerLayer + " div").css('background-color', 'transparent');
            for (index = _i = 0, _len = paths.length; _i < _len; index = ++_i) {
              path = paths[index];
              tileId = window.Map.tileIdConvert([path.x, path.y]);
              $("#tile" + Map.playerLayer + "-" + tileId).css('background', 'red');
            }
          }
          window.Map.setFocus(tileId, paths.length * 500);
          window.Map.playerTile = tileId;
          return window.Map.showDestination("#tile" + (window.Map.mapData.layers.length - 1) + "-" + (tileId - 1));
        }
      }
    };

    Map.showDestination = function(elem) {
      $('.destination').remove();
      return $(elem).append('<img class="destination" src="../resources/destination.png" />');
    };

    Map.makePath = function(toTileId) {
      var board, fromTileLoc, prop, tile, tileidc, toTileLoc, totalMapSize, x, y, _i, _j, _ref, _ref1;

      totalMapSize = this.mapData.width * this.mapData.height;
      toTileLoc = this.tileIdConvert(toTileId);
      fromTileLoc = this.tileIdConvert(Map.playerTile);
      board = [];
      for (y = _i = 0, _ref = Map.mapData.height - 1; _i < _ref; y = _i += 1) {
        board[y] = [];
        for (x = _j = 0, _ref1 = Map.mapData.width - 1; _j <= _ref1; x = _j += 1) {
          tile = this.tileIdConvert([x, y]);
          if (prop = this.tileProperties[tile]) {
            prop = prop.property;
            if (prop === 'block') {
              board[y][x] = 1;
              if (this.showBlocked) {
                tileidc = this.tileIdConvert([x, y]);
                $("#tile0-" + tileidc).css('background', 'red');
              }
            } else {
              board[y][x] = 0;
            }
          }
        }
      }
      return AStar(board, fromTileLoc, toTileLoc, 'Diagonal');
    };

    return Map;

  })();

}).call(this);
