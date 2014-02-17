var Container = require('../display/Container'),
    Sprite = require('../display/Sprite'),
    Rectangle = require('../geom/Rectangle'),
    Vector = require('../math/Vector'),
    ObjectPool = require('../utils/ObjectPool'),
    ObjectFactory = require('../utils/ObjectFactory'),
    //camera fx
    Close = require('../fx/camera/Close'),
    Fade = require('../fx/camera/Fade'),
    Flash = require('../fx/camera/Flash'),
    Scanlines = require('../fx/camera/Scanlines'),
    Shake = require('../fx/camera/Shake'),

    inherit = require('../utils/inherit'),
    math = require('../math/math'),
    C = require('../constants');

/**
 * A basic Camera object that provides some effects. It also will contain the GUI
 * to ensure they are using "screen-coords".
 *
 * @class Camera
 * @extends Container
 * @constructor
 * @param state {State} The game state this camera belongs to
 */
var Camera = function(state) {
    /**
     * The world instance this camera is tied to
     *
     * @property world
     * @type World
     */
    this.world = state.world;

    /**
     * The game instance this camera belongs to
     *
     * @property game
     * @type Game
     */
    this.game = state.game;

    /**
     * The game state this camera belongs to
     *
     * @property state
     * @type State
     */
    this.state = state;

    /**
     * The bounds of that the camera can move to
     *
     * @property bounds
     * @type Rectangle
     * @readOnly
     * @private
     */
    this.bounds = state.world.bounds.clone();

    /**
     * When following a sprite this is the space within the camera that it can move around
     * before the camera moves to track it.
     *
     * @property _deadzone
     * @type Rectangle
     * @readOnly
     * @private
     */
    this._deadzone = null;

    /**
     * The target that the camera will follow
     *
     * @property _target
     * @type Sprite
     * @readOnly
     * @private
     */
    this._target = null;

    /**
     * The target's last position, to cache if we should try to move the camera or not
     *
     * @property _targetPos
     * @type Vector
     * @readOnly
     * @private
     */
    this._targetPos = new Vector();

    /**
     * The size of the camera
     *
     * @property size
     * @type Vector
     * @readOnly
     */
    this.size = new Vector();

    /**
     * Half of the size of the camera
     *
     * @property hSize
     * @type Vector
     * @readOnly
     */
    this.hSize = new Vector();

    /**
     * The container that holds all the GUI items, direct children of Camera are effects
     *
     * @property gui
     * @type Container
     * @readOnly
     */
    this.gui = new Container();

    /**
     * An object factory for creating game objects
     *
     * @property add
     * @type ObjectFactory
     */
    this.add = new ObjectFactory(state, this.gui);

    /**
     * The fxpools for doing camera effects
     *
     * @property fxpools
     * @type Object
     * @private
     * @readOnly
     */
    this.fxpools = {
        flash: new ObjectPool(Flash, this),
        fade: new ObjectPool(Fade, this),
        shake: new ObjectPool(Shake, this),
        scanlines: new ObjectPool(Scanlines, this),
        close: new ObjectPool(Close, this)
    };

    /**
     * Flash the screen with a color. This will cover the screen in a
     * color then fade it out.
     *
     * @method flash
     * @param [color=0xFFFFFF] {Number} The color to flash the screen with
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade out
     * @param [alpha=1] {Number} The opacity of the initial flash of color (start opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Flash} The close effect that was created.
     */

    /**
     * Fade the screen into a color. This will fade into a color that will
     * eventually cover the screen.
     *
     * @method fade
     * @param [color=0xFFFFFF] {Number} The color to fade into
     * @param [duration=1000] {Number} The time it should take (in milliseconds) to fade in
     * @param [alpha=1] {Number} The opacity to fade into (final opacity)
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Fade} The close effect that was created.
     */

    /**
     * Shakes the camera around a bit.
     *
     * @method shake
     * @param [intensity=0.01] {Number} The intensity of the shaking
     * @param [duration=1000] {Number} The amount of time the screen shakes for (in milliseconds)
     * @param [direction=gf.AXIS.BOTH] {gf.AXIS} The axis to shake on
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Shake} The close effect that was created.
     */

    /**
     * Adds arcade-style scanlines to the camera viewport.
     *
     * @method scanlines - color, axis, spacing, thickness, alpha, cb
     * @param [color=0x000000] {Number} The color for the scanlines to be
     * @param [axis=gf.AXIS.HORIZONTAL] {gf.AXIS} The axis to draw the lines on
     * @param [spacing=4] {Number} Number of pixels between each line
     * @param [thickness=1] {Number} Number of pixels thick each line is
     * @param [alpha=0.3] {Number} The opacity of the lines
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Scanlines} The close effect that was created.
     */

    /**
     * Performs a "close" animation that will cover the screen with a color.
     *
     * @method close
     * @param [shape='circle'] {String} The shape to close with, can be either 'ellipse', 'circle', or 'rectangle'
     * @param [duration=1000] {Number} Number of milliseconds for the animation to complete
     * @param [position] {Vector} The position for the animation to close in on, defaults to camera center
     * @param [callback] {Function} A callback to call once the animation completes.
     * @return {fx.camera.Close} The close effect that was created.
     */

    //Dynamic addition of fx shortcuts
    var self = this;
    Object.keys(this.fxpools).forEach(function(key) {
        self[key] = function() {
            var e = self.fxpools[key].create(),
                args = Array.prototype.slice.call(arguments),
                cb = args.pop();

            if(cb !== undefined && typeof cb !== 'function')
                args.push(cb);

            args.push(this._fxCallback.bind(this, e, key, cb));

            return e.start.apply(e, args);
        };
    });

    Container.call(this);

    //add the gui child
    this.addChild(this.gui);
};

inherit(Camera, Container, {
    /**
     * The base callback for camera FX. This is called at the end of each aniamtion to
     * free the FX class back into the pool.
     *
     * @method _fxCallback
     * @param fx {mixed} The FX instance to free
     * @param type {String} The name of the instance type
     * @param [callback] {Function} The user callback to call.
     * @private
     */
    _fxCallback: function(fx, type, cb) {
        var ret;

        if(typeof cb === 'function')
            ret = cb();

        this.fxpools[type].free(fx);

        return ret;
    },
    /**
     * Follows an sprite with the camera, ensuring they are always center view. You can
     * pass a follow style to change the area an sprite can move around in before we start
     * to move with them.
     *
     * @method follow
     * @param sprite {Sprite} The sprite to follow
     * @param [style=CAMERA_FOLLOW.LOCKON] {CAMERA_FOLLOW} The style of following
     * @return {Camera} Returns itself.
     * @chainable
     */
    follow: function(spr, style) {
        if(!(spr instanceof Sprite))
            return this;

        this._target = spr;
        this._targetPos.set(null, null);

        switch(style) {
            case C.CAMERA_FOLLOW.PLATFORMER:
                var w = this.size.x / 8;
                var h = this.size.y / 3;
                this._deadzone = new Rectangle(
                    (this.size.x - w) / 2,
                    (this.size.y - h) / 2 - (h / 4),
                    w,
                    h
                );
                break;
            case C.CAMERA_FOLLOW.TOPDOWN:
                var sq4 = Math.max(this.size.x, this.size.y) / 4;
                this._deadzone = new Rectangle(
                    (this.size.x - sq4) / 2,
                    (this.size.y - sq4) / 2,
                    sq4,
                    sq4
                );
                break;
            case C.CAMERA_FOLLOW.TOPDOWN_TIGHT:
                var sq8 = Math.max(this.size.x, this.size.y) / 8;
                this._deadzone = new Rectangle(
                    (this.size.x - sq8) / 2,
                    (this.size.y - sq8) / 2,
                    sq8,
                    sq8
                );
                break;
            case C.CAMERA_FOLLOW.LOCKON:
                /* falls through */
            default:
                this._deadzone = null;
                break;
        }

        this.focusSprite(this._target);

        return this;
    },
    /**
     * Stops following any sprites
     *
     * @method unfollow
     * @return {Camera} Returns itself.
     * @chainable
     */
    unfollow: function() {
        this._target = null;
        this._targetPos.set(null, null);
        return this;
    },
    /**
     * Focuses the camera on a sprite.
     *
     * @method focusSprite
     * @param sprite {Sprite} The sprite to focus on
     * @return {Camera} Returns itself.
     * @chainable
     */
    focusSprite: function(spr) {
        var x = spr.position.x,
            y = spr.position.y,
            p = spr.parent;

        //need the transform of the sprite that doesn't take into account
        //the world object. So add up the positions not including the world position.
        while(p && p !== this.world) {
            x += p.position.x;
            y += p.position.y;
            p = p.parent;
        }

        return this.focus(
            //multiple the calculated point by the world scale for this sprite
            x * this.world.scale.x,
            y * this.world.scale.y
        );
    },
    /**
     * Focuses the camera on an x,y position. Ensures that the camera does
     * not go outside the bounds set with setBounds()
     *
     * @method focus
     * @param x {Number|Vector} The x coord to focus on, if a Vector is passed the y param is ignored
     * @param y {Number} The y coord to focus on
     * @return {Camera} Returns itself.
     * @chainable
     */
    focus: function(x, y) {
        y = x.y !== undefined ? x.y : (y || 0);
        x = x.x !== undefined ? x.x : (x || 0);

        //calculate how much we need to pan
        var goToX = x - (this.hSize.x / this.world.worldTransform.a),
            goToY = y - (this.hSize.y / this.world.worldTransform.d),
            dx = goToX + this.world.position.x, //world pos is negative
            dy = goToY + this.world.position.y;

        return this.pan(dx, dy);
    },
    /**
     * Pans the camera around by the x,y amount. Ensures that the camera does
     * not go outside the bounds set with setBounds()
     *
     * @method pan
     * @param x {Number|Vector} The x amount to pan, if a Point is passed the y param is ignored
     * @param y {Number} The y ammount to pan
     * @return {Camera} Returns itself.
     * @chainable
     */
    pan: function(dx, dy) {
        dy = dx.y !== undefined ? dx.y : (dy || 0);
        dx = dx.x !== undefined ? dx.x : (dx || 0);

        if(!dx && !dy) return;

            //world position
        var pos = this.world.position,
            //new world position
            newX = pos.x - dx,
            newY = pos.y - dy,
            b = this.bounds;

        if(b) {
            //check if X movement is illegal
            if(this._outsideBounds(-newX, -pos.y)) {
                dx = (dx < 0 ? b.x : b.right - this.size.x) + pos.x; //how far can we move since dx is too much
            }
            //check if Y movement is illegal
            if(this._outsideBounds(-pos.x, -newY)) {
                dy = (dy < 0 ? b.y : b.bottom - this.size.y) + pos.y;
            }
        }

        if(dx || dy) {
            //prevent NaN
            if(!dx) dx = 0;
            if(!dy) dy = 0;

            this.world.pan(-dx, -dy);
        }

        return this;
    },
    /**
     * Checks if a point is outside the bounds of the camera constraints.
     *
     * @method _outsideBounds
     * @param x {Number} The new X position to test
     * @param y {Number} The new Y position to test
     * @return {Boolean} true if the camera will move outside bounds to go to this point
     * @private
     */
    _outsideBounds: function(x, y) {
        //check if each corner of the camera is within the bounds
        return (
            !this.bounds.contains(x, y) || //top left
            !this.bounds.contains(x, y + this.size.y) || //bottom left
            !this.bounds.contains(x + this.size.x, y) || //top right
            !this.bounds.contains(x + this.size.x, y + this.size.y) //bottom right
        );
    },
    /**
     * Resizes the viewing area, this is called internally by your game instance
     * when you call mygame.resize(). DO NOT CALL THIS DIRECTLY
     *
     * @method resize
     * @private
     * @param w {Number} The new width
     * @param h {Number} The new height
     * @return {Camera} Returns itself.
     * @chainable
     */
    resize: function(w, h) {
        this.size.set(w, h);
        this.hSize.set(
            math.round(this.size.x / 2),
            math.round(this.size.y / 2)
        );

        return this;
    },
    /**
     * Sets the bounds the camera is allowed to go. Usually this is the world's
     * size unless you set it manually.
     *
     * @method constrain
     * @param shape {Rectangle|Polygon|Circle|Ellipse} The shape to constrain the camera into
     * @return {Camera} Returns itself.
     * @chainable
     */
    constrain: function(shape) {
        this.bounds = shape;

        return this;
    },
    /**
     * Removes the constraints of the camera, to allow free movement around the world
     *
     * @method unconstrain
     * @return {Camera} Returns itself.
     * @chainable
     */
    unconstrain: function() {
        this.bounds = null;

        return this;
    },
    /**
     * Called internally every frame. Updates all effects and the follow
     *
     * @method update
     * @param dt {Number} The delta time (in seconds) since the last update
     * @return {Camera} Returns iteself for chainability
     * @private
     */
    update: function(dt) {
        //follow sprite
        if(this._target) {
            var worldTransform = this._target.worldTransform,
                x = worldTransform.tx,
                y = worldTransform.ty;

            if(this._targetPos.x !== x || this._targetPos.y !== y) {
                this._targetPos.set(x, y);

                if(!this._deadzone) {
                    this.focusSprite(this._target);
                } else {
                    var moveX, moveY,
                        dx, dy;

                    moveX = moveY = dx = dy = 0;

                    //check less than
                    dx = x - this._deadzone.x;
                    dy = y - this._deadzone.y;

                    if(dx < 0)
                        moveX = dx;
                    if(dy < 0)
                        moveY = dy;

                    //check greater than
                    dx = x - (this._deadzone.x + this._deadzone.width);
                    dy = y - (this._deadzone.y + this._deadzone.height);

                    if(dx > 0)
                        moveX = dx;
                    if(dy > 0)
                        moveY = dy;

                    this.pan(moveX, moveY);
                }
            }
        }

        //update effects
        for(var i = 0, il = this.children.length; i < il; ++i) {
            var c = this.children[i];
            if(c.update)
                c.update(dt);
        }

        return this;
    }
});

module.exports = Camera;
