/* temporarily copied and modified from https://github.com/pfitzseb/jl-profile.js/blob/a13f2ef7852bc1782ec180e3efd4ccef2da7ba6d/dist/profile-viewer.js */
export class ProfileViewer {
    constructor(element, data, selectorLabel) {
        this.selections = [];
        this.offsetX = 0;
        this.offsetY = 0;
        this.isWheeling = false;
        this.canWheelDown = true;
        this.scrollPosition = 0;
        this.isResizing = false;
        this.isDocumentScrolling = false;
        this.isMouseMove = false;
        this.scale = window.devicePixelRatio;
        this.borderWidth = 2;
        this.padding = 2;
        this.fontConfig = '10px sans-serif';
        this.borderColor = '#fff';
        this.selectorLabel = 'Thread';
        this.boxHeight = 24;
        this.destroyed = false;
        if (typeof element === 'string') {
            element = document.querySelector(element);
        }
        if (!element) {
            throw new Error('Invalid parent element specified.');
        }
        this.container = element;
        if (selectorLabel) {
            this.selectorLabel = selectorLabel;
        }
        this.insertDOM();
        this.getStyles();
        this.registerResizeObserver();
        this.registerScrollListener();
        if (data) {
            this.setData(data);
        }
        this.getOffset();
    }
    /**
     * Remove event listeners and added child elements. The global stylesheet
     * is only removed if this is the last reference to it (i.e. there are no
     * other not-destroyed ProfileViewer instances in the DOM).
     */
    destroy() {
        this.destroyed = true;
        this.resizeObserver.disconnect();
        if (this.scrollHandler) {
            document.removeEventListener('scroll', this.scrollHandler);
        }
        if (this.stylesheet && parseInt(this.stylesheet.dataset.references) === 0) {
            document.head.removeChild(this.stylesheet);
        }
        while (this.container.firstChild) {
            this.container.removeChild(this.container.lastChild);
        }
    }
    setData(data) {
        if (this.destroyed) {
            console.error('This profile viewer is destroyed.');
            return;
        }
        if (!data) {
            this.data = data;
            this.clear();
            return;
        }
        const selections = Object.keys(data);
        selections.sort((a, b) => {
            if (a === 'all') {
                return -1;
            }
            if (b === 'all') {
                return 1;
            }
            if (a < b) {
                return -1;
            }
            if (a > b) {
                return 1;
            }
            return 0;
        });
        this.data = data;
        this.selections = selections;
        this.currentSelection = this.selections[0];
        this.activeNode = this.data[this.currentSelection];
        this.updateFilter();
        this.redraw();
    }
    setSelectorLabel(label) {
        this.selectorLabel = label;
        this.selectorLabelElement.innerText = `${label}: `;
    }
    registerCtrlClickHandler(f) {
        this.ctrlClickHandler = f;
    }
    /**
     * @deprecated Use `registerSelectionHandler` instead.
     */
    registerThreadSelectorHandler(f) {
        this.selectionHandler = f;
    }
    registerSelectionHandler(f) {
        this.selectionHandler = f;
    }
    registerScrollListener() {
        document.addEventListener('scroll', this.scrollHandler);
    }
    clear() {
        this.selections = [];
        this.currentSelection = '';
        this.activeNode = undefined;
        this.canvasCtx.clearRect(0, 0, this.canvasWidth, this.canvasHeight);
        this.hoverCanvasCtx.clearRect(0, 0, this.canvasWidth, this.canvasHeight);
    }
    isDestroyed() {
        return this.destroyed;
    }
    getStyles() {
        var _a, _b, _c;
        const style = window.getComputedStyle(this.container, null);
        const fontFamily = style.fontFamily;
        const fontSize = style.fontSize;
        this.fontConfig =
            parseInt(fontSize !== null && fontSize !== void 0 ? fontSize : '12px', 10) * this.scale +
                'px ' +
                (fontFamily !== null && fontFamily !== void 0 ? fontFamily : 'sans-serif');
        this.borderColor = (_a = style.color) !== null && _a !== void 0 ? _a : '#000';
        this.canvasCtx.font = this.fontConfig;
        this.canvasCtx.textBaseline = 'middle';
        const textMetrics = this.canvasCtx.measureText('ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]*\'"^_`abcdefghijklmnopqrstuvwxyz');
        this.boxHeight = Math.ceil((((_b = textMetrics.fontBoundingBoxDescent) !== null && _b !== void 0 ? _b : textMetrics.actualBoundingBoxDescent) +
            ((_c = textMetrics.fontBoundingBoxAscent) !== null && _c !== void 0 ? _c : textMetrics.actualBoundingBoxAscent) +
            2 * this.borderWidth +
            2 * this.padding) *
            this.scale);
        if (this.activeNode) {
            this.redraw();
        }
    }
    redraw() {
        this.canWheelDown = false;
        this.canvasCtx.clearRect(0, 0, this.canvasWidth, this.canvasHeight);
        this.clearHover();
        this.drawGraph(this.activeNode, this.canvasWidth, this.canvasHeight, 0, this.scrollPosition);
    }
    insertDOM() {
        this.insertStylesheet();
        this.canvas = document.createElement('canvas');
        this.canvas.classList.add('__profiler-canvas');
        this.canvasCtx = this.canvas.getContext('2d');
        this.hoverCanvas = document.createElement('canvas');
        this.hoverCanvas.classList.add('__profiler-hover-canvas');
        this.hoverCanvasCtx = this.hoverCanvas.getContext('2d');
        const canvasContainer = document.createElement('div');
        canvasContainer.classList.add('__profiler-canvas-container');
        canvasContainer.appendChild(this.canvas);
        canvasContainer.appendChild(this.hoverCanvas);
        canvasContainer.appendChild(this.createTooltip());
        this.container.appendChild(this.createFilterContainer());
        this.container.appendChild(canvasContainer);
        this.canvas.addEventListener('wheel', (ev) => {
            if (!this.activeNode) {
                return;
            }
            if (ev.deltaY > 0 && !this.canWheelDown) {
                return;
            }
            if (ev.deltaY < 0 && this.scrollPosition === 0) {
                if (-ev.deltaY > this.boxHeight) {
                    const parent = this.findParentNode(this.activeNode);
                    if (parent) {
                        ev.preventDefault();
                        ev.stopPropagation();
                        this.clearHover();
                        this.activeNode = parent;
                        this.redraw();
                    }
                    return;
                }
            }
            ev.preventDefault();
            ev.stopPropagation();
            if (!this.isWheeling) {
                window.requestAnimationFrame(() => {
                    this.scrollPosition = Math.min(0, this.scrollPosition - ev.deltaY);
                    this.redraw();
                    this.isWheeling = false;
                });
                this.isWheeling = true;
            }
        });
        this.canvas.addEventListener('mousemove', (ev) => {
            if (!this.isMouseMove && this.activeNode) {
                window.requestAnimationFrame(() => {
                    // XXX: this is bad
                    this.getOffset();
                    const mouseX = ev.clientX - this.offsetX;
                    const mouseY = ev.clientY - this.offsetY;
                    this.hoverCanvasCtx.clearRect(0, 0, this.canvasWidth, this.canvasHeight);
                    const didDraw = this.drawHover(this.activeNode, this.scale * mouseX, this.scale * mouseY);
                    if (didDraw) {
                        if (mouseX > this.canvasWidthCSS / 2) {
                            this.tooltipElement.style.right =
                                this.canvasWidthCSS - mouseX + 10 + 'px';
                            this.tooltipElement.style.left = 'unset';
                        }
                        else {
                            this.tooltipElement.style.right = 'unset';
                            this.tooltipElement.style.left = mouseX + 10 + 'px';
                        }
                        if (mouseY > this.canvasHeightCSS / 2) {
                            this.tooltipElement.style.bottom =
                                this.canvasHeightCSS - mouseY + 10 + 'px';
                            this.tooltipElement.style.top = 'unset';
                        }
                        else {
                            this.tooltipElement.style.bottom = 'unset';
                            this.tooltipElement.style.top = mouseY + 10 + 'px';
                        }
                        this.tooltipElement.style.display = 'block';
                    }
                    else {
                        this.tooltipElement.style.display = 'none';
                    }
                    this.isMouseMove = false;
                });
                this.isMouseMove = true;
            }
        });
        this.canvas.addEventListener('click', (ev) => {
            if (!this.activeNode) {
                return;
            }
            ev.preventDefault();
            ev.stopPropagation();
            this.getOffset();
            const mouseX = this.scale * (ev.clientX - this.offsetX);
            const mouseY = this.scale * (ev.clientY - this.offsetY);
            if (ev.ctrlKey || ev.metaKey) {
                this.runOnNodeAtMousePosition(this.activeNode, mouseX, mouseY, (node) => {
                    if (this.ctrlClickHandler) {
                        this.ctrlClickHandler(node);
                    }
                });
            }
            else {
                if (this.zoomInOnNode(this.activeNode, mouseX, mouseY)) {
                    this.scrollPosition = 0;
                    this.redraw();
                }
                else if (ev.detail === 2) {
                    // reset on double-click
                    this.resetView();
                }
            }
        });
    }
    resetView() {
        this.activeNode = this.data[this.currentSelection];
        this.scrollPosition = 0;
        this.redraw();
    }
    insertStylesheet() {
        const stylesheet = document.querySelector('#__profiler_stylesheet');
        if (stylesheet) {
            stylesheet.dataset.references = (parseInt(stylesheet.dataset.references) + 1).toString();
            this.stylesheet = stylesheet;
        }
        else {
            this.stylesheet = document.createElement('style');
            this.stylesheet.setAttribute('id', '__profiler-stylesheet');
            this.stylesheet.dataset.references = '0';
            this.stylesheet.innerText = `
                .__profiler-canvas {
                    z-index: 0;
                    position: absolute;
                    width: 100%;
                }
                .__profiler-canvas-container {
                  width: 100%;
                  height: 100%;
                  position: relative;
                }
                .__profiler-hover-canvas {
                    z-index: 1;
                    position: absolute;
                    pointer-events: none;
                    width: 100%;
                }
                .__profiler-tooltip {
                    z-index: 2;
                    display: none;
                    position: absolute;
                    background-color: #ddd;
                    border: 1px solid black;
                    padding: 5px 10px;
                    pointer-events: none;
                    max-width: 45%;
                    overflow: hidden;
                }
                .__profiler-tooltip > div {
                    line-break: anywhere;
                }
                .__profiler-filter {
                    height: 30px;
                    padding: 2px 16px;
                    margin: 0;
                    box-sizing: border-box;
                    border-bottom: 1px solid #444;
                    user-select: none;
                }
                .__profiler-reset {
                    float: right;
                }
            `;
            document.head.appendChild(this.stylesheet);
        }
    }
    createTooltip() {
        this.tooltipElement = document.createElement('div');
        this.tooltipElement.classList.add('__profiler-tooltip');
        this.tooltip = {
            count: document.createElement('span'),
            percentage: document.createElement('span'),
            function: document.createElement('code'),
            file: document.createElement('a'),
            flags: document.createElement('span'),
        };
        this.tooltip.function.classList.add('fname');
        const rows = [
            [this.tooltip.function],
            [
                this.tooltip.count,
                document.createTextNode(' ('),
                this.tooltip.percentage,
                document.createTextNode(') '),
            ],
            [this.tooltip.file],
            [this.tooltip.flags],
        ];
        for (const row of rows) {
            const rowContainer = document.createElement('div');
            for (const col of row) {
                rowContainer.appendChild(col);
            }
            this.tooltipElement.appendChild(rowContainer);
        }
        this.tooltip['ctrlClickHint'] = document.createElement('small');
        this.tooltipElement.appendChild(this.tooltip['ctrlClickHint']);
        this.container.appendChild(this.tooltipElement);
        return this.tooltipElement;
    }
    createFilterContainer() {
        this.filterContainer = document.createElement('div');
        this.filterContainer.classList.add('__profiler-filter');
        this.selectorLabelElement = document.createElement('label');
        this.selectorLabelElement.innerText = `${this.selectorLabel}: `;
        this.filterContainer.appendChild(this.selectorLabelElement);
        this.filterInput = document.createElement('select');
        this.filterInput.addEventListener('change', () => {
            this.currentSelection = this.filterInput.value;
            if (this.selectionHandler) {
                this.selectionHandler(this.currentSelection);
            }
            this.resetView();
        });
        this.filterContainer.appendChild(this.filterInput);
        const resetter = document.createElement('button');
        resetter.classList.add('__profiler-reset');
        resetter.innerText = 'reset view';
        resetter.addEventListener('click', () => {
            this.resetView();
        });
        this.filterContainer.appendChild(resetter);
        return this.filterContainer;
    }
    updateFilter() {
        while (this.filterInput.firstChild) {
            this.filterInput.removeChild(this.filterInput.lastChild);
        }
        for (const selection of this.selections) {
            const entry = document.createElement('option');
            entry.innerText = selection;
            entry.setAttribute('value', selection);
            this.filterInput.appendChild(entry);
        }
    }
    registerResizeObserver() {
        this.resizeObserver = new ResizeObserver((entries) => {
            if (!this.isResizing) {
                for (const entry of entries) {
                    if (entry.target === this.container) {
                        window.requestAnimationFrame(() => {
                            if (window.devicePixelRatio !== this.scale) {
                                this.scale = window.devicePixelRatio;
                                this.getStyles();
                            }
                            this.canvasWidth = Math.round(entry.contentRect.width * this.scale);
                            this.canvasHeight = Math.round((entry.contentRect.height - 30) * this.scale);
                            this.canvasWidthCSS = entry.contentRect.width;
                            this.canvasHeightCSS = entry.contentRect.height;
                            this.canvas.width = this.canvasWidth;
                            this.canvas.height = this.canvasHeight;
                            this.hoverCanvas.width = this.canvasWidth;
                            this.hoverCanvas.height = this.canvasHeight;
                            this.redraw();
                            this.isResizing = false;
                        });
                    }
                }
                this.isResizing = true;
            }
        });
        this.resizeObserver.observe(this.container);
    }
    scrollHandler(e) {
        if (!this.isDocumentScrolling) {
            window.requestAnimationFrame(() => {
                this.getOffset();
                this.isDocumentScrolling = false;
            });
            this.isDocumentScrolling = true;
        }
    }
    getOffset() {
        const box = this.canvas.getBoundingClientRect();
        this.offsetX = box.left;
        this.offsetY = box.top;
    }
    // hash of function named, used to seed PRNG
    nodeHash(node) {
        const hashString = node.file + node.line;
        let hash = 0;
        for (let i = 0; i < hashString.length; i++) {
            const char = hashString.charCodeAt(i);
            hash = (hash << 5) - hash + char;
            hash = hash & hash;
        }
        return hash;
    }
    // Simple PRNG from https://stackoverflow.com/a/47593316/12113178
    mulberry32(a) {
        return function () {
            let t = (a += 0x6d2b79f5);
            t = Math.imul(t ^ (t >>> 15), t | 1);
            t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
            return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
        };
    }
    // modifies the normal color by three stable random values drawn from a
    // PRNG seeded by the node hash
    modifyNodeColorByHash(r, g, b, hash, range = 255) {
        const rng = this.mulberry32(hash);
        if (r === g && g === b) {
            r = g = b = Math.min(255, Math.max(0, r + (rng() - 0.5) * range));
        }
        else {
            r = Math.min(range, Math.max(0, r + (rng() - 0.5) * range));
            g = Math.min(range, Math.max(0, g + (rng() - 0.5) * range));
            b = Math.min(range, Math.max(0, b + (rng() - 0.5) * range));
        }
        return {
            r,
            g,
            b,
        };
    }
    modifyNodeColorByCount(r, g, b, count, range = 255) {
            r = Math.min(255, g * (1-count) );
            g = Math.min(255, g * (1-count) );
            b = Math.min(255, b * (1-count) );
        return {
            r,
            g,
            b,
        };
    }
    nodeColors(node, hash) {
        let r, g, b;
        let a = 1;
        /**
        if (node.flags & 0x01) {
            // runtime-dispatch
            ;
            ({ r, g, b } = this.modifyNodeColorByHash(204, 103, 103, hash, 20));
        }
        else if (node.flags & 0x02) {
            // gc
            ;
            ({ r, g, b } = this.modifyNodeColorByHash(204, 153, 68, hash, 20));
        }
        else if (node.flags & 0x08) {
            // compilation?
            ;
            ({ r, g, b } = this.modifyNodeColorByHash(100, 100, 100, hash, 60));
        }
        else {
            // default
            ;
            ({ r, g, b } = this.modifyNodeColorByHash(64, 99, 221, hash));
        }
        if (node.flags & 0x10) {
            // C frame
            a = 0.5;
        }
        */

        // if (node.count > 0) {
        // better performance
        ;
        if (node.count_change / node.count > 100) {
            ({ r, g, b } = this.modifyNodeColorByCount(0, 0, 0, node.count_change / node.count ));
        }
        else if (node.count_change/ node.count > 0) {
            ({ r, g, b } = this.modifyNodeColorByCount(255, 255, 255, node.count_change / node.count ));
            r = 255 - 20;
            g = g - 20;
            b = b - 20;

        }
        else {
            ({ r, g, b } = this.modifyNodeColorByCount(255, 255, 255, - node.count_change / node.count ));
            b = 255 - 20;
            g = g - 20;
            r = r - 20;
        }

        return {
            fill: 'rgba(' + r + ',' + g + ',' + b + ',' + a + ')',
            stroke: 'rgba(' + 0.8 * r + ',' + 0.8 * g + ',' + 0.8 * b + ',' + a + ')',
            text: 'rgba(255, 255, 255, ' + Math.max(0.6, a) + ')',
        };
    }
    drawGraph(node, width, height, x, y) {
        if (!node) {
            return;
        }
        this.canvasCtx.font = this.fontConfig;
        this.canvasCtx.textBaseline = 'middle';
        if (y + this.boxHeight >= 0) {
            const hash = this.nodeHash(node);
            const { fill, stroke, text } = this.nodeColors(node, hash);
            this.drawNode(node.func, fill, stroke, text, width, x, y);
        }
        node.pos = {
            x,
            y,
            width,
            height: this.boxHeight,
        };
        if (y + this.boxHeight <= this.canvasHeight) {
            for (const child of node.children) {
                const w = width * (child.fraction || child.count / node.count);
                this.drawGraph(child, w, height, x, y + this.boxHeight);
                x += w;
            }
        }
        else {
            this.canWheelDown = true;
        }
    }
    drawNode(text, color, bColor, textColor, width, x, y) {
        if (width < 1) {
            width = 1;
        }
        const drawBorder = false; //width > 20*this.borderWidth;
        this.canvasCtx.fillStyle = color;
        this.canvasCtx.beginPath();
        this.canvasCtx.rect(x, y + this.borderWidth, width, this.boxHeight - this.borderWidth);
        this.canvasCtx.closePath();
        this.canvasCtx.fill();
        if (drawBorder) {
            this.canvasCtx.fillStyle = bColor;
            this.canvasCtx.beginPath();
            this.canvasCtx.rect(x, y + this.borderWidth, this.borderWidth, this.boxHeight - this.borderWidth);
            this.canvasCtx.closePath();
            this.canvasCtx.fill();
        }
        const textWidth = width - 2 * this.padding - 2 * this.borderWidth;
        if (textWidth > 10) {
            this.canvasCtx.save();
            this.canvasCtx.beginPath();
            this.canvasCtx.rect(x + this.borderWidth + this.padding, y + this.borderWidth + this.padding, textWidth, this.boxHeight - this.borderWidth - 2 * this.padding);
            this.canvasCtx.closePath();
            this.canvasCtx.clip();
            this.canvasCtx.fillStyle = textColor;
            this.canvasCtx.fillText(text, x + this.borderWidth + this.padding, y + this.boxHeight / 2 + this.borderWidth);
            this.canvasCtx.restore();
        }
    }
    updateTooltip(node) {
        this.tooltip.function.innerText = node.func;
        if (node.file || node.line > 0) {
            this.tooltip.file.innerText = node.file + ':' + node.line;
        }
        else {
            this.tooltip.file.innerText = '';
        }
        this.tooltip.count.innerText = (node.countLabel || (node.count + ' samples')).toString();
        let percentageText = ((100 * node.count) /
            this.data[this.currentSelection].count).toFixed() + '% of root';
        if (this.activeNode.count != this.data[this.currentSelection].count) {
            percentageText = percentageText + ', ' + ((100 * node.count) /
                this.activeNode.count).toFixed() + '% of selection';
        }
        this.tooltip.percentage.innerText = percentageText;
        const flags = [];
        if (node.flags & 0x01) {
            flags.push('runtime-dispatch');
        }
        if (node.flags & 0x02) {
            flags.push('GC');
        }
        if (node.flags & 0x08) {
            flags.push('compilation');
        }
        if (node.flags & 0x10) {
            flags.push('task');
        }
        let flagString = '';
        if (flags.length > 0) {
            flagString = 'Flags: ' + flags.join(', ');
        }
        this.tooltip.flags.innerText = flagString;
        if (this.ctrlClickHandler) {
            this.tooltip['ctrlClickHint'].innerText =
                'Ctrl/Cmd+Click to open this file';
        }
    }
    drawHoverNode(node) {
        this.hoverCanvasCtx.fillStyle = this.borderColor;
        this.hoverCanvasCtx.fillRect(node.pos.x, node.pos.y + this.borderWidth, Math.max(1, node.pos.width), node.pos.height - this.borderWidth);
        const innerWidth = node.pos.width - this.borderWidth * 2 * this.scale;
        if (innerWidth > 1) {
            this.hoverCanvasCtx.clearRect(node.pos.x + this.borderWidth * this.scale, node.pos.y + 2 * this.borderWidth * this.scale, innerWidth, node.pos.height - this.borderWidth * 3 * this.scale);
        }
        this.updateTooltip(node);
    }
    clearHover() {
        this.hoverCanvasCtx.clearRect(0, 0, this.canvasWidth, this.canvasHeight);
        this.tooltipElement.style.display = 'none';
    }
    drawHover(node, mouseX, mouseY) {
        let found = false;
        this.runOnNodeAtMousePosition(node, mouseX, mouseY, (node) => {
            this.drawHoverNode(node);
            found = true;
        });
        return found;
    }
    runOnNodeAtMousePosition(root, x, y, f) {
        if (x >= Math.floor(root.pos.x) &&
            x <= Math.ceil(root.pos.x + root.pos.width) &&
            y >= root.pos.y) {
            if (y <= root.pos.y + root.pos.height) {
                f(root);
                return true;
            }
            else {
                for (const child of root.children) {
                    if (this.runOnNodeAtMousePosition(child, x, y, f)) {
                        return true;
                    }
                }
            }
        }
        return false;
    }
    zoomInOnNode(node, mouseX, mouseY) {
        let found = false;
        this.runOnNodeAtMousePosition(node, mouseX, mouseY, (node) => {
            this.clearHover();
            this.activeNode = node;
            found = true;
        });
        return found;
    }
    // ideally this wouldn't require tree traversal at all
    findParentNode(target, current = null) {
        if (current === null) {
            current = this.data[this.currentSelection];
        }
        for (const child of current.children) {
            if (child === target) {
                return current;
            }
            else {
                const found = this.findParentNode(target, child);
                if (found) {
                    return found;
                }
            }
        }
        return null;
    }
}
//# sourceMappingURL=profile-viewer.js.map