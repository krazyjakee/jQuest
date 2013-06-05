(function() {
  this.keyframes = (function() {
    function keyframes() {}

    keyframes.frameCollection = [];

    keyframes.cache = [0];

    keyframes.expando = 'data' + +new Date();

    keyframes.data = function(elem) {
      var cacheIndex, nextCacheIndex;

      cacheIndex = elem[this.expando];
      nextCacheIndex = this.cache.length;
      if (!cacheIndex) {
        cacheIndex = elem[this.expando] = nextCacheIndex;
        this.cache[cacheIndex] = {};
      }
      return {
        get: function(key) {
          return keyframes.cache[cacheIndex][key];
        },
        set: function(key, val) {
          keyframes.cache[cacheIndex][key] = val;
          return val;
        }
      };
    };

    keyframes.styleElem = '';

    keyframes.browserCode = '';

    keyframes.setBrowserCode = function() {
      var ua;

      ua = navigator.userAgent;
      if (ua.indexOf('Opera') !== -1) {
        return '-o-';
      } else if (ua.indexOf('MSIE') !== -1) {
        return '-ms-';
      } else if (ua.indexOf('WebKit') !== -1) {
        return '-webkit-';
      } else if (navigator.product === 'Gecko') {
        return '-moz-';
      } else {
        return '';
      }
    };

    keyframes.generate = function() {
      var browserType, css, frameData, frameName, key, list, _ref, _ref1;

      keyframes.styleElem.innerHTML = '';
      browserType = this.browserCode;
      list = [];
      _ref = this.frameCollection;
      for (key in _ref) {
        frameName = _ref[key];
        css = '@' + browserType + 'keyframes ' + frameName.data.name + '{';
        _ref1 = frameName.data;
        for (key in _ref1) {
          frameData = _ref1[key];
          if (key !== 'name') {
            css += key + '{';
            css += frameData + ';}';
          }
        }
        css += '}\n';
        list.push(css);
      }
      list.push(' .boostKeyframe{transform:scale3d(1,1,1);}');
      return this.styleElem.innerHTML = list.join('');
    };

    keyframes.removeHead = function() {
      return this.styleElem.innerHTML = '';
    };

    keyframes.reset = function(elem, callback) {
      var animationkey;

      animationkey = this.browserCode + 'animation';
      elem.style[animationkey + '-play-state'] = 'running';
      elem.style[animationkey] = 'none';
      this.data(elem).set("keyframe", false);
      clearInterval(this.data(elem).get('keyframeTimer'));
      clearTimeout(this.data(elem).get('keyframeTimer'));
      if (callback) {
        return setTimeout(callback, 0);
      }
    };

    keyframes.pause = function(elem) {
      elem.style[this.browserCode + 'animation-play-state'] = 'paused';
      clearInterval(this.data(elem).get('keyframeTimer'));
      return clearTimeout(this.data(elem).get('keyframeTimer'));
    };

    keyframes.resume = function(elem) {
      return elem.style[this.browserCode + 'animation-play-state'] = 'running';
    };

    keyframes.add = function(frameData) {
      var data, kfname, _i, _len, _results;

      _results = [];
      for (_i = 0, _len = frameData.length; _i < _len; _i++) {
        data = frameData[_i];
        kfname = data.name;
        _results.push(this.frameCollection[kfname] = {
          data: data
        });
      }
      return _results;
    };

    keyframes.playNextItem = function(elem, callback) {
      var frames;

      frames = keyframes.data(elem).get('keyframeQueue');
      if (frames.length) {
        keyframes.play(elem, frames[0], function() {
          return keyframes.playNextItem(elem, callback);
        });
        return keyframes.data(elem).set('keyframeQueue', frames.slice(1));
      } else {
        if (callback) {
          return callback(elem);
        }
      }
    };

    keyframes.play = function(elem, frameOptions, callback, reset, regen) {
      var _playKeyframe;

      if (reset == null) {
        reset = true;
      }
      if (regen == null) {
        regen = true;
      }
      _playKeyframe = function() {
        var animationcss, animationkey, delay, duration, frameOptSplit, key, name, opt, repeat;

        if (typeof frameOptions === 'string') {
          frameOptions = frameOptions.trim();
          frameOptSplit = frameOptions.split(' ');
          name = frameOptSplit[0];
          duration = parseInt(frameOptSplit[1]);
          delay = parseInt(frameOptSplit[3]);
          repeat = parseInt(frameOptSplit[4]);
          frameOptSplit[1] += 'ms';
          frameOptSplit[3] += 'ms';
          animationcss = frameOptSplit.join(' ');
        } else if (frameOptions.length) {
          keyframes.data(elem).set('keyframeQueue', frameOptions);
          keyframes.playNextItem(elem, callback);
          return;
        } else {
          name = frameOptions.name;
          duration = frameOptions.duration;
          delay = frameOptions.delay;
          repeat = frameOptions.repeat;
          frameOptions.duration += 'ms';
          frameOptions.delay += 'ms';
          animationcss = '';
          for (key in frameOptions) {
            opt = frameOptions[key];
            animationcss += opt + ' ';
          }
          animationcss = animationcss.trim();
        }
        animationkey = keyframes.browserCode + 'animation';
        if (regen) {
          keyframes.generate();
        }
        if (repeat !== 'infinite') {
          if (callback) {
            keyframes.data(elem).set('keyframeTimer', setTimeout(function() {
              return callback(elem);
            }, (duration + delay) * repeat));
          }
          setTimeout(function() {
            return keyframes.data(elem).set('keyframe', false);
          }, (duration + delay) * repeat);
        } else {
          if (callback) {
            keyframes.data(elem).set('keyframeTimer', setTimeout(function() {
              callback(elem);
              return keyframes.data(elem).set('keyframeTimer', setInterval(function() {
                return callback(elem);
              }, duration));
            }, duration + delay));
          }
        }
        elem.setAttribute("style", keyframes.browserCode + 'animation-play-state: running; ' + animationkey + ': ' + animationcss);
        return keyframes.data(elem).set('keyframe', name);
      };
      if (reset) {
        return this.reset(elem, _playKeyframe);
      } else {
        return _playKeyframe();
      }
    };

    return keyframes;

  })();

  document.getElementsByTagName('head')[0].innerHTML += '<style id="keyframes-style" type="text/css"></style>';

  keyframes.styleElem = document.getElementById('keyframes-style');

  keyframes.browserCode = keyframes.setBrowserCode();

}).call(this);
