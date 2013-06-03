(function() {
  window.Character = (function() {
    function Character() {}

    Character.settings = {
      fileExtension: 'json',
      spriteDirectory: ''
    };

    Character.sprites = [];

    Character.characters = [];

    Character.playerSprite = false;

    Character.playerMoving = false;

    Character.loadSprite = function(charSource, callback) {
      var spritefile;

      spritefile = "" + (this.settings.spriteDirectory + charSource) + "." + this.settings.fileExtension;
      return $.getJSON(spritefile, function(json) {
        window.Character.sprites[json.id] = json;
        if (callback) {
          return callback(json.id);
        }
      });
    };

    Character.placeCharacter = function(id, sprite, location) {
      var charid, loc;

      loc = window.Map.tileIdPosition(location);
      $("#layer" + Map.playerLayer).append("<div class=\"character\" id=\"character" + id + "\" />");
      charid = "#character" + id;
      return $(charid).width(this.sprites[sprite].cellWidth).height(this.sprites[sprite].cellHeight).css("background-image", "url(" + this.sprites[sprite].sprite + ")").css("background-position", "-" + this.sprites[sprite].cellWidth + "px 0px").css("left", loc[0]).css("top", loc[1]);
    };

    Character.loadPlayer = function(sprite) {
      window.Character.playerSprite = sprite;
      window.Character.characters[0] = {
        sprite: sprite,
        timer: []
      };
      return window.Character.placeCharacter(0, sprite, window.Map.playerTile);
    };

    Character.playerMove = function(paths) {
      /*
          If the character is already moving at the end location is redefined, it's best to remove the first
          of the next path as animation becomes more fluid.
      */

      var index, _i, _ref;

      if (this.playerMoving) {
        paths = paths.slice(1);
        this.stopAnimation(0);
        for (index = _i = 0, _ref = this.characters[0].timer.length; 0 <= _ref ? _i < _ref : _i > _ref; index = 0 <= _ref ? ++_i : --_i) {
          clearTimeout(this.characters[0].timer[index]);
        }
      }
      $('#character0').clearQueue().stop();
      this.playerMoving = true;
      if (paths[0]) {
        this.animate(0, paths[0].direction);
      }
      return $.each(paths, function(index, path) {
        var loc;

        loc = Map.tileIdPosition([path.x, path.y]);
        if (path.direction) {
          window.Character.characters[0].timer[index] = setTimeout(function() {
            return window.Character.animate(0, path.direction);
          }, index * 310);
        }
        return $('#character0').animate({
          left: loc[0],
          top: loc[1]
        }, 300, 'linear', function() {
          if (index === paths.length - 1) {
            window.Character.playerMoving = false;
            window.Character.stopAnimation(0);
            $('.destination').remove();
          }
          return Map.playerTile = Map.tileIdConvert([path.x + 1, path.y]);
        });
      });
    };

    Character.animate = function(id, direction) {
      window.Character.stopAnimation(id);
      return window.Character.characters[id].animation = {
        timer: setInterval(function() {
          var x, y;

          if (window.Character.characters[id].animation.step > 5) {
            window.Character.characters[id].animation.step = 0;
          }
          y = 0;
          switch (direction) {
            case "n":
              y = 3;
              break;
            case "e":
              y = 1;
              break;
            case "s":
              y = 0;
              break;
            case "w":
              y = 2;
              break;
            case "ne":
              y = 6;
              break;
            case "nw":
              y = 7;
              break;
            case "se":
              y = 5;
              break;
            case "sw":
              y = 4;
          }
          y *= window.Character.sprites[id].cellHeight;
          x = window.Character.characters[id].animation.step * window.Character.sprites[id].cellWidth;
          $("#character" + id).css("background-position", "-" + x + "px -" + y + "px");
          return window.Character.characters[id].animation.step++;
        }, 200),
        step: 0,
        y: 0
      };
    };

    Character.stopAnimation = function(id) {
      var bgpos;

      if (window.Character.characters[id].animation) {
        clearInterval(window.Character.characters[id].animation.timer);
        bgpos = $("#character" + id).css('background-position').split(' ');
        bgpos[0] = "-" + this.sprites[id].cellWidth + "px";
        return $("#character" + id).css('background-position', bgpos.join(' '));
      }
    };

    return Character;

  })();

}).call(this);
